"""Figure 2 追加パネル: マッピング率の検体内対応比較（参照 → FVBパーソナルゲノム）"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

BASE = "/usr/local/jupyter/FVB_B6_glo"
OUT = f"{BASE}/results/figures"
C_B6, C_FVB = "#8172B2", "#C44E52"

t = pd.read_csv(f"{BASE}/results/tables/Table2_sequencing_mapping.csv")
REF = "Uniquely mapped, reference (%)"
CB = "Uniquely mapped, cond. B (%)"

fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2),
                          gridspec_kw=dict(width_ratios=[1.25, 1]))

# --- 左: before-after ---
ax = axes[0]
for _, r in t.iterrows():
    c = C_FVB if r.Strain == "FVB" else C_B6
    ax.plot([0, 1], [r[REF], r[CB]], "-o", color=c, ms=5, lw=1.2, alpha=.8,
            markeredgecolor="white", markeredgewidth=.6)
for st, c in [("B6", C_B6), ("FVB", C_FVB)]:
    d = t[t.Strain == st]
    ax.plot([0, 1], [d[REF].mean(), d[CB].mean()], "-", color=c, lw=3.2, alpha=.95,
            solid_capstyle="round", zorder=5)
ax.set_xlim(-0.28, 1.42); ax.set_xticks([0, 1])
ax.set_xticklabels(["GRCm39\n(reference = B6)", "FVB personalised\ngenome (cond. B)"], fontsize=11.9)
ax.set_ylabel("Uniquely mapped reads (%)", fontsize=12.6)
ax.set_title("A  Mapping rate, within-sample comparison", fontsize=14, weight="bold")
ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=11.2)
for st, c, dy in [("FVB", C_FVB, 0.10), ("B6", C_B6, -0.14)]:
    d = t[t.Strain == st]
    ax.text(1.06, d[CB].mean() + dy, st, color=c, fontsize=13.3, weight="bold", va="center")
ax.annotate("FVBP1-3", xy=(0, 88.28), xytext=(-0.24, 89.6), fontsize=9.5, color="#666",
            arrowprops=dict(arrowstyle="-", color="#999", lw=.7))

# --- 右: Δ の分布 ---
ax = axes[1]
for i, (st, c) in enumerate([("B6", C_B6), ("FVB", C_FVB)]):
    d = t[t.Strain == st]
    delta = (d[CB] - d[REF]).values
    ax.scatter(np.random.normal(i, .055, len(delta)), delta, s=42, color=c,
               alpha=.85, edgecolor="white", linewidth=.7, zorder=3)
    ax.plot([i - .26, i + .26], [delta.mean()] * 2, color="black", lw=2, zorder=4)
    p = stats.ttest_rel(d[CB], d[REF]).pvalue
    ax.text(i, delta.mean() + (0.085 if delta.mean() > 0 else -0.16),
            f"{delta.mean():+.2f} pt\n$p$={p:.0e}", ha="center", fontsize=10.6)
ax.axhline(0, color="#888", lw=.9, ls="--")
ax.set_xlim(-.55, 1.55); ax.set_xticks([0, 1]); ax.set_xticklabels(["B6", "FVB"], fontsize=12.6)
ax.set_ylabel("Change in mapping rate (pt)", fontsize=12.6)
ax.set_ylim(-0.75, 0.75)
ax.set_title("B  Strain-specific effect", fontsize=14, weight="bold")
ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=11.2)
ax.text(.5, .50, "all 6 FVB samples improve;\nall 6 B6 samples worsen",
        transform=ax.transAxes, ha="center", va="center", fontsize=10.4, color="#444",
        bbox=dict(boxstyle="round,pad=.35", fc="#f7f7f7", ec="#ccc", lw=.5))

fig.tight_layout()
fig.savefig(f"{OUT}/Figure2_panelAB_mapping.png", dpi=600)
fig.savefig(f"{OUT}/Figure2_panelAB_mapping.pdf")
print(f"保存: {OUT}/Figure2_panelAB_mapping.png / .pdf")
