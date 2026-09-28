#!/usr/bin/env python3
"""Figure 3C–F (bias correction). Same layout and styling as the 'Figure 2' block of plot_omics.py
(the file that an earlier version of the analysis produced as Figure3_panelsCF), but every number is
read from the bias-QC outputs instead of being written into the script:
  --metrics  bias_metrics.tsv written by rerun2/A/scripts/bias_qc_param.py
  --qc-dir   directory with log2FC_{E13.5,P1}.tsv from the same run (scatter in panel E)
  --matched  Figure3F_matched_shifts.tsv, from rerun5/density (panel F)
The panel-D note reports Spearman rho for both conditions; the former 'indel-driven' note is removed.

v2 (2026-09-28): the fourth panel is replaced. It used to plot the reverse-control shift against
1.99, the mean |log2 fold change| between FVB/N and C57BL/6J in genes chosen for having a large
strain difference. Those are not the same quantity, so that comparison is withdrawn. The panel now
shows the forward and reverse shifts as the same quantity on the same genes, which is what the
manuscript reports."""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--metrics", required=True); ap.add_argument("--qc-dir", required=True)
ap.add_argument("--matched", required=True, help="source_data/Figure3F_matched_shifts.tsv")
ap.add_argument("--out-prefix", required=True); ap.add_argument("--letters", default="ABCD")
a = ap.parse_args()
M = pd.read_csv(a.metrics, sep="\t", index_col=0)
m = lambda k, tp: float(M.loc[k, tp])
L = a.letters
C_FVB = "#C44E52"; COL = {"E13.5": "#4C72B0", "P1": "#55A868"}
def style(ax): ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=11.2)

fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.6))

# A: convergence in high-bias genes
ax = axes[0, 0]
x = np.arange(3); w = 0.36
for i, tp in enumerate(["E13.5", "P1"]):
    v = [m(f"high_bias_mean_abs_{arm}", tp) for arm in ["uncorrected", "condA", "condB"]]
    ax.bar(x + (i - 0.5) * w, v, w, label=tp, color=COL[tp], alpha=.85)
ax.set_xticks(x); ax.set_xticklabels(["Uncorrected", "Cond. A\n(SNP only)", "Cond. B\n(SNP+indel)"], fontsize=12.9)
ax.set_ylabel("Mean |log2FC|", fontsize=14.5); ax.set_ylim(1.5, 2.1)
ax.set_title(f"{L[0]}  Convergence in high-bias genes", fontsize=16.1, weight="bold")
ax.legend(fontsize=12.9, frameon=False); style(ax)
ax.text(.5, .03, f"E13.5 n={int(m('high_bias_n','E13.5'))}, P1 n={int(m('high_bias_n','P1'))} genes (|uncorrected log2FC|>1)",
        transform=ax.transAxes, ha="center", fontsize=10.9, color="#555")

# B: correction magnitude vs variant density quartile
ax = axes[0, 1]
Q = ["Q1_low", "Q2", "Q3", "Q4_high"]; xq = np.arange(4)
vals = {}
for tp, mk in [("E13.5", "o"), ("P1", "s")]:
    for arm, ls, al, lw, ms in [("condB", "-", 1.0, 1.6, 5), ("condA", "--", .45, 1.2, 4)]:
        y = [m(f"quartile_{q}_{arm}", tp) for q in Q]; vals[(tp, arm)] = y
        ax.plot(xq, y, mk + ls, color=COL[tp], alpha=al, lw=lw, ms=ms, label=f"{tp} cond. {arm[-1]}")
ymax = max(max(v) for v in vals.values())
ax.set_xticks(xq); ax.set_xticklabels(["Q1\n(low)", "Q2", "Q3", "Q4\n(high)"], fontsize=12.9)
ax.set_xlabel("FVB variant density quartile", fontsize=14.5); ax.set_ylabel("Mean |correction|", fontsize=14.5)
ax.set_title(f"{L[1]}  Correction magnitude vs variant density", fontsize=16.1, weight="bold")
ax.set_ylim(0, ymax * 1.45)
h, lab = ax.get_legend_handles_labels(); order = [0, 2, 1, 3]
ax.legend([h[i] for i in order], [lab[i] for i in order], fontsize=11.3, frameon=False, loc="upper center", ncol=2); style(ax)
ax.text(0, -0.34, "Spearman ρ (|correction| vs density)\n"
        f"cond. A: {m('spearman_rho_condA','E13.5'):.3f} (E13.5), {m('spearman_rho_condA','P1'):.3f} (P1)\n"
        f"cond. B: {m('spearman_rho_condB','E13.5'):.3f} (E13.5), {m('spearman_rho_condB','P1'):.3f} (P1)",
        transform=ax.transAxes, ha="left", va="top", fontsize=10.6, color="#444",
        bbox=dict(boxstyle="round,pad=.28", fc="#fdf6e3", ec="#ddd", lw=.5))

# C: condition A vs B agreement
ax = axes[1, 0]
for tp in ["E13.5", "P1"]:
    d = pd.read_csv(os.path.join(a.qc_dir, f"log2FC_{tp}.tsv"), sep="\t", index_col=0)
    sub = d[["log2FC_condA", "log2FC_condB"]].dropna()
    sub = sub.sample(min(3000, len(sub)), random_state=1)
    ax.scatter(sub.iloc[:, 0], sub.iloc[:, 1], s=3, alpha=.18, color=COL[tp], label=tp, rasterized=True)
ax.plot([-4, 4], [-4, 4], "k--", lw=.8, alpha=.5)
ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
ax.set_xlabel("log2FC (condition A)", fontsize=14.5); ax.set_ylabel("log2FC (condition B)", fontsize=14.5)
ax.set_title(f"{L[2]}  Condition A vs B agreement", fontsize=16.1, weight="bold")
ax.legend(fontsize=12.9, frameon=False, markerscale=3, loc="lower right"); style(ax)
ax.text(.03, .95, f"Pearson r\nE13.5 {m('r_condA_condB','E13.5'):.2f}\nP1 {m('r_condA_condB','P1'):.2f}", transform=ax.transAxes, va="top", fontsize=11.3,
        bbox=dict(boxstyle="round,pad=.3", fc="#f5f5f5", ec="#bbb", lw=.5))

# D (published panel F): matched forward and reverse shifts, same quantity on the same genes
ax = axes[1, 1]
MS = pd.read_csv(a.matched, sep="\t")
shown = MS[(MS.timepoint == "E13.5") & (MS.arm == "condA")].iloc[0]
bars = [("Forward (FVB/N)", float(shown.forward), C_FVB), ("Reverse (C57BL/6J)", float(shown.reverse), "#8172B2")]
for i, (lab_, v, col) in enumerate(bars):
    ax.bar(i, v, width=.5, color=col, alpha=.9)
    ax.text(i, v + 0.0008, f"{v:.3f}", ha="center", fontsize=12.1)
ax.set_xticks(range(len(bars))); ax.set_xticklabels([b[0] for b in bars], fontsize=12.1)
ax.set_xlim(-0.7, 1.7); ax.set_ylim(0, max(b[1] for b in bars) * 1.35)
ax.set_ylabel("Mean absolute log$_2$ shift", fontsize=14.5)
ax.set_title(f"{L[3]}  Matched forward and reverse shifts", fontsize=16.1, weight="bold"); style(ax)
ax.text(0, 1.10, f"{shown.timepoint}, condition {shown.arm[-1]}; same {int(shown.n_genes):,} genes",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=11.3, color="#555")
ratio = MS.forward / MS.reverse
rho_lo, rho_hi = MS.spearman_signed.max(), MS.spearman_signed.min()
nv = MS[(MS.timepoint == "E13.5") & (MS.arm == "condA")].iloc[0]
ax.text(0, -0.30, "Across both stages and genome conditions:\n"
        f"Forward / reverse magnitude: approximately {ratio.min():.1f}\u2013{ratio.max():.1f}\n"
        f"Signed shifts: Spearman \u03c1 = {rho_lo:.2f} to {rho_hi:.2f}\n"
        f"No FVB variant: forward {nv.novariant_forward:.3f}; reverse {nv.novariant_reverse:.3f}",
        transform=ax.transAxes, ha="left", va="top", fontsize=10.6, color="#555")

fig.tight_layout(); fig.subplots_adjust(hspace=0.78)
fig.savefig(a.out_prefix + ".png", dpi=600, bbox_inches="tight"); fig.savefig(a.out_prefix + ".pdf", bbox_inches="tight"); plt.close(fig)
vals_df = pd.DataFrame({"panel": ["A"] * 6 + ["B"] * 16 + ["B"] * 4 + ["C"] * 2,
    "item": [f"{tp} {arm}" for tp in ["E13.5", "P1"] for arm in ["uncorrected", "condA", "condB"]]
            + [f"{tp} {arm} {q}" for tp in ["E13.5", "P1"] for arm in ["condB", "condA"] for q in Q]
            + [f"rho {arm} {tp}" for arm in ["condA", "condB"] for tp in ["E13.5", "P1"]]
            + ["r E13.5", "r P1"],
    "value": [m(f"high_bias_mean_abs_{arm}", tp) for tp in ["E13.5", "P1"] for arm in ["uncorrected", "condA", "condB"]]
            + [v for tp in ["E13.5", "P1"] for arm in ["condB", "condA"] for v in vals[(tp, arm)]]
            + [m(f"spearman_rho_{arm}", tp) for arm in ["condA", "condB"] for tp in ["E13.5", "P1"]]
            + [m("r_condA_condB", "E13.5"), m("r_condA_condB", "P1")]})
d_rows = []
for r in MS.itertuples():
    d_rows += [{"panel": "D", "item": f"{r.timepoint} {r.arm} forward", "value": r.forward},
               {"panel": "D", "item": f"{r.timepoint} {r.arm} reverse", "value": r.reverse},
               {"panel": "D", "item": f"{r.timepoint} {r.arm} forward/reverse", "value": round(r.forward / r.reverse, 4)},
               {"panel": "D", "item": f"{r.timepoint} {r.arm} Spearman rho (signed)", "value": round(r.spearman_signed, 4)}]
d_rows += [{"panel": "D", "item": "no-variant genes, E13.5 cond. A forward", "value": round(float(nv.novariant_forward), 4)},
           {"panel": "D", "item": "no-variant genes, E13.5 cond. A reverse", "value": round(float(nv.novariant_reverse), 4)}]
pd.concat([vals_df, pd.DataFrame(d_rows)], ignore_index=True).to_csv(a.out_prefix + "_plotted_values.tsv", sep="\t", index=False)
print("saved", a.out_prefix)
