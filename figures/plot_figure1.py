"""表現型データの図（Figure 1候補）"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tost_options import tost as tost_opt
ap = argparse.ArgumentParser()
ap.add_argument("--data", required=True); ap.add_argument("--outdir", required=True)
ap.add_argument("--tost-mode", choices=["pooled", "welch"], required=True)
ap.add_argument("--seed", type=int, default=1, help="seed for the horizontal jitter (the original script had none)")
args = ap.parse_args(); OUT = args.outdir; os.makedirs(OUT, exist_ok=True); np.random.seed(args.seed)

df = pd.read_csv(args.data)
df["glom_per_g_kw"] = df.glomeruli_per_kidney / df.kidney_weight_g
df["glom_per_g_bw"] = df.glomeruli_per_kidney / df.body_weight_g
df["kw_bw_pct"] = df.kidney_weight_g / df.body_weight_g * 100

ORDER = ["BALB", "DBA", "B6", "FVB"]
FULL = {"BALB": "BALB/c", "DBA": "DBA/2", "B6": "C57BL/6J", "FVB": "FVB/N"}
COL = {"BALB": "#4C72B0", "DBA": "#55A868", "B6": "#8172B2", "FVB": "#C44E52"}

def holm_bonferroni(pvals):
    """Holm–Bonferroni 補正（多重比較）。"""
    pvals = np.asarray(pvals)
    order = np.argsort(pvals)
    n = len(pvals)
    adj = np.empty(n)
    running_max = 0.0
    for rank, i in enumerate(order):
        running_max = max(running_max, (n - rank) * pvals[i])
        adj[i] = min(running_max, 1.0)
    return adj

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tukey import tukey_hsd

STAR_LOG = []

def stars(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"

def panel(ax, col, ylabel, title, fmt="{:.0f}"):
    for i, s in enumerate(ORDER):
        d = df[df.strain == s]
        for sx, mk in [("M", "o"), ("F", "^")]:
            v = d[d.sex == sx][col]
            ax.scatter(np.random.normal(i, 0.07, len(v)), v, s=26, marker=mk,
                       color=COL[s], alpha=.75, edgecolor="white", linewidth=.6, zorder=3)
        m, sd = d[col].mean(), d[col].std()
        ax.plot([i - .28, i + .28], [m, m], color="black", lw=1.8, zorder=4)
        ax.errorbar(i, m, yerr=sd, color="black", capsize=4, lw=1.1, zorder=4)
    ax.set_xticks(range(4)); ax.set_xticklabels([FULL[s] for s in ORDER], fontsize=11.2)
    ax.set_ylabel(ylabel, fontsize=12.6); ax.set_title(title, fontsize=14, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11.2)

    # 有意差マーク: Tukey HSD（4系統すべての対比較）から FVB/N 対 他3系統を取る。
    # 以前は Student t 検定 + Holm 補正で、本文が引用する Tukey HSD と別の検定だった。
    # 絶対糸球体数の FVB/N 対 C57BL/6J で結論が食い違い（t+Holm で **、Tukey では ns）、
    # 本文の「絶対数の差は多重比較補正後に有意でない」という記述と矛盾していた。
    others = ["B6", "DBA", "BALB"]  # FVBに近い列(x=2)から遠い列(x=0)の順＝ブラケットが交差しない順
    tk = tukey_hsd(df[col].values, df.strain.values)
    padj = [next(v for k, v in tk.items() if set(k) == {"FVB", s}) for s in others]
    for s, p in zip(others, padj):
        STAR_LOG.append(dict(panel=title.split()[0], measure=col, comparison=f"FVB vs {s}",
                             test="Tukey HSD (all four strains)", tukey_p=p, star=stars(p)))
    y0, y1 = ax.get_ylim()
    step = (y1 - y0) * 0.115
    top = y1 + step * 0.4
    for k, (s, p) in enumerate(zip(others, padj)):
        xi = ORDER.index(s)
        y = top + k * step
        ax.plot([xi, xi, 3, 3], [y, y + step * 0.18, y + step * 0.18, y],
                lw=1.1, color="black", clip_on=False)
        ax.text((xi + 3) / 2, y + step * 0.18, stars(p), ha="center", va="bottom",
                fontsize=11.2, clip_on=False)
    ax.set_ylim(y0, top + len(others) * step + step * 0.6)

fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.6))

panel(axes[0, 0], "glomeruli_per_kidney", "Glomeruli per kidney", "A  Absolute glomerular number")
panel(axes[0, 1], "glom_per_g_kw", "Glomeruli / g kidney weight", "B  Normalised to kidney mass")
panel(axes[1, 0], "kw_bw_pct", "Kidney weight / body weight (%)", "C  Relative kidney mass")

# D: 内部対照 — FVBとBALBは腎重量が同等だが糸球体数は約2倍違う
ax = axes[1, 1]
for s in ORDER:
    d = df[df.strain == s]
    a = .55 if s in ("DBA", "B6") else .95
    for sx, mk in [("M", "o"), ("F", "^")]:
        dd = d[d.sex == sx]
        ax.scatter(dd.kidney_weight_g, dd.glomeruli_per_kidney, s=34, marker=mk, color=COL[s],
                   alpha=a, edgecolor="white", linewidth=.6,
                   label=FULL[s] if sx == "M" else None, zorder=3)
ax.set_xlabel("Kidney weight (g)", fontsize=12.6)
ax.set_ylabel("Glomeruli per kidney", fontsize=12.6)
ax.set_title("D  Internal control: FVB vs BALB/c", fontsize=14, weight="bold")
ax.legend(fontsize=10.5, frameon=False, loc="upper left")
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(labelsize=11.2)

# 内部対照の統計: 「差がない」は TOST 同等性検定で示す（t検定のp>0.05では主張できない）
m = df[df.sex == "M"]

kw = tost_opt(m[m.strain == "FVB"].kidney_weight_g, m[m.strain == "BALB"].kidney_weight_g, args.tost_mode)
kb = tost_opt(m[m.strain == "FVB"].kw_bw_pct, m[m.strain == "BALB"].kw_bw_pct, args.tost_mode)
g_p = stats.ttest_ind(m[m.strain == "FVB"].glomeruli_per_kidney,
                       m[m.strain == "BALB"].glomeruli_per_kidney).pvalue
dflab = "pooled TOST, df = 8" if args.tost_mode == "pooled" else "Welch TOST"
kbp = "< 0.001" if kb["TOST_P"] < 0.001 else f"={kb['TOST_P']:.3f}"
ax.text(0, -0.34,
        f"FVB/N vs BALB/c (males; {dflab})\n"
        f"kidney weight   Δ={kw['diff']:+.3f} g, 90% CI [{kw['CI90_lo']:+.3f}, {kw['CI90_hi']:+.3f}]\n"
        f"                          margin ±{kw['margin']:.3f} g, TOST $p$={kw['TOST_P']:.3f}\n"
        f"KW/BW              Δ={kb['diff']:+.3f} pp, 90% CI [{kb['CI90_lo']:+.3f}, {kb['CI90_hi']:+.3f}]\n"
        f"                          margin ±{kb['margin']:.3f} pp, TOST $p$ {kbp}\n"
        f"glomeruli          $p$={g_p:.3f}  (differ ~2-fold)\n"
        f"pp, percentage points; margin ±20% of BALB/c mean",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.8, linespacing=1.35,
        bbox=dict(boxstyle="round,pad=.35", fc="#f5f5f5", ec="#bbbbbb", lw=.6))
pd.DataFrame([dict(measure="kidney_weight", **kw), dict(measure="kw_bw_ratio", **kb)]).to_csv(f"{OUT}/Figure1_plotted_values_TOST.tsv", sep="\t", index=False)

fig.text(.005, .5, "circles = males, triangles = females; bars = mean ± SD (n=5 per sex)",
         rotation=90, va="center", fontsize=9.8, color="#555555")
fig.tight_layout(rect=[.02, 0, 1, 1])
fig.savefig(f"{OUT}/Figure1.png", dpi=600, bbox_inches="tight")
fig.savefig(f"{OUT}/Figure1.pdf", bbox_inches="tight")
import csv as _csv
with open(f"{OUT}/figure1_stars_tukey.tsv", "w", newline="") as _f:
    _w = _csv.DictWriter(_f, fieldnames=list(STAR_LOG[0]), delimiter="\t")
    _w.writeheader(); _w.writerows(STAR_LOG)
print("saved", args.tost_mode, "-", len(STAR_LOG), "star comparisons ->", f"{OUT}/figure1_stars_tukey.tsv")
