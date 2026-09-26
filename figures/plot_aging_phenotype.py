"""加齢コホートの表現型データ（Figure 2）

★ 週齢ラベル（10wk / 8wk 等）と n 数は Methods 確定後に検証すること。
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

BASE = "/usr/local/jupyter/FVB_B6_glo/submission"
OUT = f"{BASE}/figures"

df = pd.read_csv(f"{BASE}/source_data/aging_cohort_raw.csv")

FULL = {"B6": "C57BL/6J", "FVB": "FVB/N"}
COL = {"B6": "#8172B2", "FVB": "#C44E52"}


def two_way_anova(sub, log=False):
    """strain x age の二元配置分散分析。全項の表を返す。"""
    d = sub.copy()
    d["y"] = np.log(d.value) if log else d.value
    return sm.stats.anova_lm(smf.ols("y ~ C(strain) * C(age_label)", data=d).fit(), typ=2)


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


def stars(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"


def trajectory(ax, measure, ages, ylabel, title, log_interaction=False, age_labels=None):
    """系統ごとの加齢推移。個体点 + 平均±SD を折れ線で結ぶ。
    age_labels: 表示用のラベル（省略時は ages をそのまま表示）。
    データ抽出は常に ages（実際の age_label 値）で行う。
    各週齢で B6 vs FVB の有意差マークを付け、統計サマリーは軸の下に表示する。"""
    sub = df[df.measure == measure]
    vals = {"B6": {}, "FVB": {}}
    for s in ["B6", "FVB"]:
        means, sds, xs = [], [], []
        for i, a in enumerate(ages):
            v = sub[(sub.strain == s) & (sub.age_label == a)].value
            vals[s][i] = v
            ax.scatter(np.random.normal(i, 0.055, len(v)), v, s=22, color=COL[s],
                       alpha=.55, edgecolor="white", linewidth=.5, zorder=3)
            means.append(v.mean()); sds.append(v.std()); xs.append(i)
        ax.errorbar(xs, means, yerr=sds, color=COL[s], lw=1.8, marker="o", ms=5,
                    capsize=3.5, zorder=4, label=FULL[s])
    ax.set_xticks(range(len(ages))); ax.set_xticklabels(age_labels or ages, fontsize=11.2)
    ax.set_ylabel(ylabel, fontsize=12.6)
    ax.set_title(title, fontsize=14, weight="bold", loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11.2)
    ax.legend(fontsize=10.5, frameon=False)

    # 各週齢での B6 vs FVB（Welch t検定 + Holm補正）を軸上部に * で表示
    raw_p = []
    for i in range(len(ages)):
        b6, fvb = vals["B6"][i], vals["FVB"][i]
        if log_interaction:
            b6, fvb = np.log(b6), np.log(fvb)
        raw_p.append(stats.ttest_ind(fvb, b6, equal_var=False).pvalue)
    padj = holm_bonferroni(raw_p)
    y0, y1 = ax.get_ylim()
    pad = (y1 - y0) * 0.09
    for i in range(len(ages)):
        top_i = max(vals["B6"][i].max(), vals["FVB"][i].max())
        ax.text(i, top_i + pad * 0.2, stars(padj[i]), ha="center", va="bottom",
                fontsize=11.9, weight="bold", clip_on=False)
    ax.set_ylim(y0, y1 + pad)

    # 統計サマリー（系統の主効果・週齢の主効果・交互作用）を軸の下に表示
    aov = two_way_anova(sub, log=log_interaction)
    p_s = aov.loc["C(strain)", "PR(>F)"]
    p_a = aov.loc["C(age_label)", "PR(>F)"]
    p_i = aov.loc["C(strain):C(age_label)", "PR(>F)"]
    suffix = " (log)" if log_interaction else ""
    ax.text(0, -0.30,
            f"two-way ANOVA{suffix}\n"
            f"strain $p$ = {p_s:.2g}\n"
            f"age $p$ = {p_a:.2g}\n"
            f"strain × age $p$ = {p_i:.2g}",
            transform=ax.transAxes, ha="left", va="top", fontsize=9.5, color="#333333",
            linespacing=1.4,
            bbox=dict(boxstyle="round,pad=.3", fc="#f5f5f5", ec="#bbbbbb", lw=.5))
    return p_i


def two_group(ax, measure, ylabel, title, box_xy=(.97, .04), box_ha="right", box_va="bottom"):
    """単一時点での2系統比較。t検定と Mann-Whitney を併記。現在は未使用。"""
    sub = df[df.measure == measure]
    vals = {}
    for i, s in enumerate(["B6", "FVB"]):
        v = sub[sub.strain == s].value
        vals[s] = v
        ax.scatter(np.random.normal(i, 0.06, len(v)), v, s=28, color=COL[s],
                   alpha=.8, edgecolor="white", linewidth=.6, zorder=3)
        ax.plot([i - .25, i + .25], [v.mean()] * 2, color="black", lw=1.8, zorder=4)
        ax.errorbar(i, v.mean(), yerr=v.std(), color="black", capsize=4, lw=1.1, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels([FULL["B6"], FULL["FVB"]], fontsize=11.2)
    ax.set_ylabel(ylabel, fontsize=12.6)
    ax.set_title(title, fontsize=14, weight="bold", loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11.2)
    pt = stats.ttest_ind(vals["FVB"], vals["B6"], equal_var=False).pvalue
    pu = stats.mannwhitneyu(vals["FVB"], vals["B6"]).pvalue
    # 境界域であることを図上で明示する（t検定のみを示さない）
    ax.text(*box_xy, f"Welch $p$ = {pt:.3f}\nMann–Whitney $p$ = {pu:.3f}",
            transform=ax.transAxes, ha=box_ha, va=box_va, fontsize=9.1,
            linespacing=1.3, color="#333333",
            bbox=dict(boxstyle="round,pad=.3", fc="#f5f5f5", ec="#bbbbbb", lw=.5))


# ★ 糸球体数（絶対数・腎重量あたり密度）は年齢間で比較していない。
#   篩分画法の回収率が腎サイズに依存し、8週→1年で絶対数が約20%「増加」する
#   （ネフロン数は生後3-4日で固定されるため生物学的にありえない）。
#   生データは source_data/aging_cohort_raw.csv に収録してあるが図示しない。

CROP_DIR = f"{OUT}/source_images/cropped"
CROP_UM = 150.0          # prep_pas_crops.py と一致させること
BAR_UM = 50.0            # 図中に描き直すスケールバー

fig = plt.figure(figsize=(13.2, 9.3))
gs = fig.add_gridspec(2, 12, height_ratios=[1.0, 1.06], hspace=.26, wspace=2.2)

# --- A: 代表PAS像（2系統 × 2週齢）---
panels = [("B6_PAS_1", "C57BL/6J", "10 wk"), ("B6_PAS_2", "C57BL/6J", "1 year"),
          ("FVB_PAS_1", "FVB/N", "10 wk"), ("FVB_PAS_2", "FVB/N", "1 year")]
for i, (fname, strain, age) in enumerate(panels):
    ax = fig.add_subplot(gs[0, i * 3:(i + 1) * 3])
    ax.imshow(mpimg.imread(f"{CROP_DIR}/{fname}_crop.png"),
              extent=[0, CROP_UM, 0, CROP_UM])
    ax.set_xticks([]); ax.set_yticks([])
    c = COL["B6"] if strain == "C57BL/6J" else COL["FVB"]
    for sp in ax.spines.values():
        sp.set_edgecolor(c); sp.set_linewidth(1.6)
    ax.set_title(f"{strain}   {age}", fontsize=12.6, color=c, weight="bold", pad=4)
    # スケールバーを描き直す（元画像の焼き込みはクロップで除外済み）
    x0, y0 = CROP_UM - BAR_UM - 8, 8
    ax.add_patch(Rectangle((x0, y0), BAR_UM, 2.6, color="black", zorder=5))
    if i == 0:
        ax.text(x0 + BAR_UM / 2, y0 + 4.5, f"{BAR_UM:.0f} µm", ha="center",
                va="bottom", fontsize=10.5, zorder=5)

fig.text(.013, .975, "A  Representative glomeruli (PAS)", fontsize=14.7,
         weight="bold", va="top", ha="left")

# --- B–D: 定量 ---
axB = fig.add_subplot(gs[1, 0:4])
axC = fig.add_subplot(gs[1, 4:8])
axD = fig.add_subplot(gs[1, 8:12])

trajectory(axB, "kidney_weight", ["8wk", "1yr", "1.5yr"],
           "Kidney weight (g)", "B  Kidney mass with age")
trajectory(axC, "glom_diameter", ["10wk", "1yr"],
           "Glomerular long-axis diameter (µm)", "C  Glomerular hypertrophy")
trajectory(axD, "uacr", ["10wk", "52wk"],
           "UACR (mg/g creatinine)", "D  Urinary albumin excretion", log_interaction=True,
           age_labels=["10wk", "1yr"])

fig.text(.055, .012,
         "Male mice only. B–D: points, individual animals; bars, mean ± SD. "
         "A: all four images at identical magnification.",
         ha="left", va="bottom", fontsize=10.1, color="#555555")
fig.subplots_adjust(left=.055, right=.99, top=.89, bottom=.23)
fig.savefig(f"{OUT}/Figure2.png", dpi=600, bbox_inches="tight")
fig.savefig(f"{OUT}/Figure2.pdf", bbox_inches="tight")
print(f"保存: {OUT}/Figure2.png / .pdf")

# ★ ポドサイト数・WT1（8週齢, n=5/群）は今回の投稿では報告しない。
#   元データは input/論文用データ_7_03.xlsx に保持されている。
#   再度含める場合は source_data に行を戻し、two_group() で作図すること。
