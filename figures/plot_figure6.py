#!/usr/bin/env python3
"""Figure 6: candidate prioritisation and the evaluation of two additional strategies.

v3 (2026-09-28): rebuilt as the four-panel figure of the submitted manuscript. The previous
version drew three panels (motif-disruption forest, window sensitivity, concordance) and
recomputed each one from the raw intermediates with bedtools. It predated two corrections:
the funnel was not shown at all, and its concordance panel used the superseded background
definition. This version draws all four panels from the aggregated tables in source_data/,
so it runs without the analysis directory and plots exactly the numbers that were reported:

  A  funnel            1,236 -> 92 -> 7 -> 5   (source_data/Figure6A_funnel.tsv,
                                                from rerun5/ties/funnel_comparison.tsv,
                                                row "all ties kept")
  B  motif disruption  SIX2 0.78, OSR1 4.10 (P = 0.0037), WT1 0.75, over 1,592 private SNVs
                       in H3K27ac regions    (source_data/Figure6B_odds_ratios.tsv)
  C  window size       SIX2 odds ratio at +/-100 to +/-2,000 bp
                       (source_data/Figure6C_window_sensitivity.tsv)
  D  concordance       48,401 / 48,510 (99.8%) and 15 / 15 (100%)
                       (source_data/Figure6D_concordance_5genes.tsv)

Panels B and C are the values the earlier script computed (rerun3/figures/F7/panel{A,B}_values.tsv);
only its stale caption, which quoted the pre-allele-wise variant count, is corrected here.

Usage: plot_figure6.py --funnel F --odds F --windows F --concordance F --n-variants N --out-prefix P
"""
import argparse
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

ap = argparse.ArgumentParser()
ap.add_argument("--funnel", required=True, help="source_data/Figure6A_funnel.tsv")
ap.add_argument("--odds", required=True, help="source_data/Figure6B_odds_ratios.tsv")
ap.add_argument("--windows", required=True, help="source_data/Figure6C_window_sensitivity.tsv")
ap.add_argument("--concordance", required=True, help="source_data/Figure6D_concordance_5genes.tsv")
ap.add_argument("--n-variants", type=int, default=1592,
                help="private SNVs in H3K27ac regions that panel B stratifies (rows of rerun3/F7/allreg_full.tsv)")
ap.add_argument("--out-prefix", required=True)
a = ap.parse_args()

C_BLUE, C_RED, C_GREY = "#4C72B0", "#C44E52", "#999999"
def style(ax):
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=10.5)

fig, axes = plt.subplots(2, 2, figsize=(13.6, 11.0),
                         gridspec_kw={"wspace": 0.30, "hspace": 0.42})

# ---------------- A: prioritisation funnel ----------------
ax = axes[0][0]; ax.set_axis_off()
fn = pd.read_csv(a.funnel, sep="\t").sort_values("step")
shades = ["#DCE4F0", "#BACBE4", "#7D9DC9", "#2F5597"]
top, step_h, gap = 1.0, 0.145, 0.075
cx = 0.24
half = [0.22, 0.187, 0.154, 0.121]
for i, (r, col) in enumerate(zip(fn.itertuples(), shades)):
    y1 = top - i * (step_h + gap); y0 = y1 - step_h
    w0, w1 = half[i], half[i] * 0.86
    ax.add_patch(Polygon([(cx - w0, y1), (cx + w0, y1), (cx + w1, y0), (cx - w1, y0)],
                         closed=True, facecolor=col, edgecolor="none",
                         transform=ax.transAxes, clip_on=False))
    ax.text(cx, (y0 + y1) / 2, f"{int(r.n_genes):,}", transform=ax.transAxes,
            ha="center", va="center", fontsize=16.5, weight="bold",
            color="white" if i >= 2 else "#1F3864")
    ax.text(0.50, (y0 + y1) / 2, str(r.label).replace("\\n", "\n"),
            transform=ax.transAxes, ha="left", va="center", fontsize=9.8)
    if i < len(fn) - 1:
        ax.annotate("", xy=(cx, y0 - gap + 0.012), xytext=(cx, y0 - 0.008),
                    xycoords=ax.transAxes, textcoords=ax.transAxes,
                    arrowprops=dict(arrowstyle="-|>", color="#777777", lw=1.1))
genes = str(fn.candidates_display_order.dropna().iloc[0]).split(";")
y_last = top - (len(fn) - 1) * (step_h + gap) - step_h
ax.text(0.0, y_last - 0.085, "Final genes: " + ", ".join(genes), transform=ax.transAxes,
        ha="left", va="top", fontsize=10.8, weight="bold")
ax.text(0.0, y_last - 0.145, "Gene-level counts; all nearest-gene ties retained.",
        transform=ax.transAxes, ha="left", va="top", fontsize=9.2, color="#666666")
ax.set_title("A  Candidate-gene prioritisation", fontsize=13.2, weight="bold", loc="left", x=-0.02)

# ---------------- B: motif-disruption enrichment ----------------
ax = axes[0][1]; style(ax)
od = pd.read_csv(a.odds, sep="\t").set_index("TF").loc[["SIX2", "OSR1", "WT1"]].reset_index()
for i, r in od.iterrows():
    ax.plot([r.CI_lo, r.CI_hi], [i, i], color="#444444", lw=1.6)
    ax.plot([r.OR], [i], marker="s", ms=9, color=C_BLUE, zorder=3)
    p = f"{r.P:.4f}".rstrip("0") if r.P < 0.01 else f"{r.P:.2f}"
    ax.text(1.06, (i + 0.5) / len(od), f"OR = {r.OR:.2f}\np = {p}",
            transform=ax.transAxes, ha="left", va="center", fontsize=10.2)
ax.axvline(1, color=C_RED, ls="--", lw=1.2)
ax.set_yticks(range(len(od))); ax.set_yticklabels(od.TF, fontsize=12.0)
ax.set_xscale("log"); ax.set_xlim(0.22, 20)
ax.set_xticks([0.25, 0.5, 1, 2, 4, 8, 16])
ax.set_xticklabels(["0.25", "0.5", "1", "2", "4", "8", "16"], fontsize=10.5)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_ylim(-0.7, len(od) - 0.3)
ax.set_xlabel("Odds ratio (95% CI)", fontsize=12.0)
ax.set_title("B  Motif-disruption enrichment", fontsize=13.2, weight="bold", loc="left", x=-0.06)
ax.text(0.0, -0.20, f"Variants in H3K27ac regions (n = {a.n_variants:,});\n"
                    "stratified by TF-specific motif-disruption calls.",
        transform=ax.transAxes, ha="left", va="top", fontsize=9.4, color="#777777")

# ---------------- C: window-size sensitivity ----------------
ax = axes[1][0]; style(ax)
wn = pd.read_csv(a.windows, sep="\t").sort_values("window_bp")
x = range(len(wn))
ax.plot(x, wn.OR, "o-", color=C_BLUE, lw=1.9, ms=7.5)
ax.axhline(1, color=C_RED, ls="--", lw=1.2)
ax.set_xticks(list(x)); ax.set_xticklabels([f"±{int(w)}" for w in wn.window_bp], fontsize=10.5)
ax.set_xlabel("Window around ChIP-seq summit (bp)", fontsize=12.0)
ax.set_ylabel("Odds ratio (SIX2)", fontsize=12.0)
ax.set_ylim(0.38, 1.62)
ax.set_yticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6])
ax.set_title("C  Insensitive to window size", fontsize=13.2, weight="bold", loc="left", x=-0.10)
for i, r in enumerate(wn.itertuples()):
    below = i == 0
    ax.annotate(f"p = {r.P:.2f}", (i, r.OR), textcoords="offset points",
                xytext=(0, -20 if below else 13), ha="center", fontsize=9.6, color="#555555")

# ---------------- D: four-strain concordance ----------------
ax = axes[1][1]; style(ax)
cd = pd.read_csv(a.concordance, sep="\t")
bars = [(f"Private variants with\ncomplete genotypes\n(n = {int(cd.n.iloc[0]):,})",
         float(cd.pct.iloc[0]), f"{int(cd.concordant.iloc[0]):,} / {int(cd.n.iloc[0]):,}", C_GREY),
        (f"Prioritised regulatory\nvariants\n(n = {int(cd.n.iloc[1]):,})",
         float(cd.pct.iloc[1]), f"{int(cd.concordant.iloc[1])} / {int(cd.n.iloc[1])}", C_BLUE)]
for i, (lab, pct, inner, col) in enumerate(bars):
    ax.bar(i, pct, color=col, alpha=.9, width=.5)
    ax.text(i, pct + 2.0, f"{pct:.0f}%" if float(pct).is_integer() else f"{pct:.1f}%",
            ha="center", fontsize=12.2, weight="bold")
    ax.text(i, pct * 0.62, inner, ha="center", fontsize=9.2,
            color="#555555" if col == C_GREY else "white")
ax.set_xticks(range(len(bars))); ax.set_xticklabels([b[0] for b in bars], fontsize=9.6)
ax.set_xlim(-0.7, 1.7); ax.set_ylim(0, 112)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_ylabel("Concordance with kidney-mass-\nnormalised phenotype rank (%)", fontsize=10.6)
ax.set_title("D  Four-strain concordance", fontsize=13.2, weight="bold", loc="left", x=-0.14)
ax.text(0.5, -0.30, "This criterion is non-discriminatory.", transform=ax.transAxes,
        ha="center", va="top", fontsize=9.2, color="#777777")

for ext in ("png", "pdf"):
    fig.savefig(f"{a.out_prefix}.{ext}", dpi=300, bbox_inches="tight")

vals = ([{"panel": "A", "item": str(r.label).replace("\\n", " "), "value": int(r.n_genes)} for r in fn.itertuples()]
        + [{"panel": "A", "item": "final genes", "value": ", ".join(genes)}]
        + [{"panel": "B", "item": f"{r.TF} OR / P", "value": f"{r.OR:.4f} / {r.P:.4g}"} for _, r in od.iterrows()]
        + [{"panel": "B", "item": "variants stratified", "value": a.n_variants}]
        + [{"panel": "C", "item": f"SIX2 OR at +/-{int(r.window_bp)} bp", "value": f"{r.OR:.4f} (P = {r.P:.4f})"} for r in wn.itertuples()]
        + [{"panel": "D", "item": str(cd.bar.iloc[i]), "value": f"{int(cd.concordant.iloc[i])}/{int(cd.n.iloc[i])} = {cd.pct.iloc[i]}%"} for i in range(len(cd))])
pd.DataFrame(vals).to_csv(f"{a.out_prefix}_plotted_values.tsv", sep="\t", index=False)
print(f"saved {a.out_prefix}")
