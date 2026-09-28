#!/usr/bin/env python3
"""Figure 3A-B: within-sample mapping rate on GRCm39 versus on the FVB/NJ pseudo-genome.

v2 (2026-09-28): the table is an argument and defaults to the final requantification
(source_data/Table2_sequencing_mapping_final.csv) rather than a fixed path to the superseded
results/tables/Table2_sequencing_mapping.csv. Both give the same panel, because the mapping-rate
columns are unchanged, but only one of them is the table the manuscript reports. Panel labels
follow the submitted figure.

Usage: plot_figure3_panelsAB.py --table F --out-prefix P
"""
import argparse
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--table", required=True, help="source_data/Table2_sequencing_mapping_final.csv")
ap.add_argument("--out-prefix", required=True)
a = ap.parse_args()
C_B6, C_FVB = "#8172B2", "#C44E52"

t = pd.read_csv(a.table)
REF = "Uniquely mapped, reference (%)"
CB = "Uniquely mapped, cond. B (%)"

fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2),
                          gridspec_kw=dict(width_ratios=[1.25, 1]))

# --- left: before / after, one line per sample ---
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
ax.set_xticklabels(["GRCm39\n(reference)", "FVB/NJ pseudo-genome\n(condition B)"], fontsize=11.9)
ax.set_ylabel("Uniquely mapped reads (%)", fontsize=12.6)
ax.set_title("A  Within-sample mapping rates", fontsize=14, weight="bold")
ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=11.2)
for st, c, dy, lab in [("B6", C_B6, 0.22, "C57BL/6J"), ("FVB", C_FVB, -0.22, "FVB/N")]:
    d = t[t.Strain == st]
    ax.text(1.06, d[CB].mean() + dy, lab, color=c, fontsize=12.2, weight="bold", va="center")
ax.annotate("FVBP1-3", xy=(0, 88.28), xytext=(-0.24, 89.6), fontsize=9.5, color="#666",
            arrowprops=dict(arrowstyle="-", color="#999", lw=.7))

VALS = []
# --- right: distribution of the within-sample change ---
ax = axes[1]
for i, (st, c) in enumerate([("B6", C_B6), ("FVB", C_FVB)]):
    d = t[t.Strain == st]
    delta = (d[CB] - d[REF]).values
    ax.scatter(np.random.normal(i, .055, len(delta)), delta, s=42, color=c,
               alpha=.85, edgecolor="white", linewidth=.7, zorder=3)
    ax.plot([i - .26, i + .26], [delta.mean()] * 2, color="black", lw=2, zorder=4)
    p = stats.ttest_rel(d[CB], d[REF]).pvalue
    ax.text(i, delta.mean() + (0.085 if delta.mean() > 0 else -0.16),
            f"{delta.mean():+.3f} points\n$P$ = {p / 10 ** int(np.floor(np.log10(p))):.1f} \u00d7 10$^{{{int(np.floor(np.log10(p)))}}}$",
            ha="center", fontsize=10.2)
    VALS.append(dict(panel="B", strain={"B6": "C57BL/6J", "FVB": "FVB/N"}[st],
                     mean_change_points=round(float(delta.mean()), 4), paired_t_P=float(p), n=len(delta)))
ax.axhline(0, color="#888", lw=.9, ls="--")
ax.set_xlim(-.55, 1.55); ax.set_xticks([0, 1]); ax.set_xticklabels(["C57BL/6J", "FVB/N"], fontsize=12.6)
ax.set_ylabel("Change in mapping rate (percentage points)", fontsize=12.1)
ax.set_ylim(-0.75, 0.75)
ax.set_title("B  Strain-specific change", fontsize=14, weight="bold")
ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=11.2)
ax.text(.5, .50, "All six FVB/N samples improve;\nall six C57BL/6J samples worsen.",
        transform=ax.transAxes, ha="center", va="center", fontsize=10.4, color="#444",
        bbox=dict(boxstyle="round,pad=.35", fc="#f7f7f7", ec="#ccc", lw=.5))

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{a.out_prefix}.{ext}", dpi=300, bbox_inches="tight")
pd.DataFrame(VALS).to_csv(f"{a.out_prefix}_plotted_values.tsv", sep="\t", index=False)
print(f"saved {a.out_prefix}")
