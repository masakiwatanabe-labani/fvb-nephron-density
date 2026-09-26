#!/usr/bin/env python3
"""Figure 3C–F (bias correction). Same layout and styling as the 'Figure 2' block of plot_omics.py
(the file submitted as Figure3_panelsCF), but every number is read from the bias-QC outputs instead of
being written into the script:
  --metrics  bias_metrics.tsv written by rerun2/A/scripts/bias_qc_param.py
  --qc-dir   directory with log2FC_{E13.5,P1}.tsv from the same run (scatter in panel C)
The panel-B note reports Spearman rho for both conditions; the former 'indel-driven' note is removed."""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--metrics", required=True); ap.add_argument("--qc-dir", required=True)
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
ax.text(0, -0.24, "Spearman ρ (|correction| vs density)\n"
        f"cond. A: {m('spearman_rho_condA','E13.5'):.2f} (E13.5), {m('spearman_rho_condA','P1'):.2f} (P1)\n"
        f"cond. B: {m('spearman_rho_condB','E13.5'):.2f} (E13.5), {m('spearman_rho_condB','P1'):.2f} (P1)",
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

# D: reverse control
ax = axes[1, 1]
lab = ["E13.5\ncond. A", "E13.5\ncond. B", "P1\ncond. A", "P1\ncond. B"]
val = [m("reverse_mean_abs_condA", "E13.5"), m("reverse_mean_abs_condB", "E13.5"), m("reverse_mean_abs_condA", "P1"), m("reverse_mean_abs_condB", "P1")]
for i, (v, tp, al) in enumerate(zip(val, ["E13.5", "E13.5", "P1", "P1"], [.5, .9, .5, .9])):
    ax.bar(i, v, color=COL[tp], alpha=al)
fwd = m("high_bias_mean_abs_uncorrected", "E13.5")
ax.axhline(fwd, color=C_FVB, ls="--", lw=1.2)
ax.set_xticks(range(4)); ax.set_xticklabels(lab, fontsize=12.1)
ax.set_ylabel("Mean |log2FC|", fontsize=14.5); ax.set_yscale("log"); ax.set_ylim(0.01, 5)
ax.set_title(f"{L[3]}  Reverse check: B6 reads → FVB genome", fontsize=16.1, weight="bold"); style(ax)
ax.text(.5, fwd * 1.3, f"forward bias in high-bias genes ({fwd:.2f})", ha="center", fontsize=11.3, color=C_FVB)

fig.tight_layout(); fig.subplots_adjust(hspace=0.62)
fig.savefig(a.out_prefix + ".png", dpi=600, bbox_inches="tight"); fig.savefig(a.out_prefix + ".pdf", bbox_inches="tight"); plt.close(fig)
pd.DataFrame({"panel": ["A"] * 6 + ["B"] * 16 + ["B"] * 4 + ["C"] * 2 + ["D"] * 5,
    "item": [f"{tp} {arm}" for tp in ["E13.5", "P1"] for arm in ["uncorrected", "condA", "condB"]]
            + [f"{tp} {arm} {q}" for tp in ["E13.5", "P1"] for arm in ["condB", "condA"] for q in Q]
            + [f"rho {arm} {tp}" for arm in ["condA", "condB"] for tp in ["E13.5", "P1"]]
            + ["r E13.5", "r P1"] + ["E13.5 condA", "E13.5 condB", "P1 condA", "P1 condB", "forward E13.5"],
    "value": [m(f"high_bias_mean_abs_{arm}", tp) for tp in ["E13.5", "P1"] for arm in ["uncorrected", "condA", "condB"]]
            + [v for tp in ["E13.5", "P1"] for arm in ["condB", "condA"] for v in vals[(tp, arm)]]
            + [m(f"spearman_rho_{arm}", tp) for arm in ["condA", "condB"] for tp in ["E13.5", "P1"]]
            + [m("r_condA_condB", "E13.5"), m("r_condA_condB", "P1")] + val + [fwd]}).to_csv(a.out_prefix + "_plotted_values.tsv", sep="\t", index=False)
print("saved", a.out_prefix)
