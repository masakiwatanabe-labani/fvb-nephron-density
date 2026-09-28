#!/usr/bin/env python3
"""Figure 2: age-associated kidney and glomerular phenotypes in C57BL/6J and FVB/N.

v2 (2026-09-28): the tests are brought into line with the manuscript. The previous version used one
test for all three quantitative panels (Welch's t-test with Holm correction across the ages within a
panel, on log-transformed values for UACR). Three things changed:

  B  kidney weight   Tukey-Kramer over the six strain-by-age groups, which is the family the
                     manuscript reports, rather than Welch plus Holm within the panel.
  C  glomerular      unchanged: Welch's t-test, Holm-corrected across the two ages.
     diameter
  D  UACR            between-strain comparisons are two-sided Welch's t-tests on the UNTRANSFORMED
                     values and are not corrected, because they are the values the manuscript
                     reports. The log transform is used only for the two-way ANOVA. The former
                     log-plus-Holm comparison is removed. Mann-Whitney U is reported as a two-sided
                     exact sensitivity analysis; the one-sided values that the earlier draft quoted
                     (0.0096 and 0.0027) are not used, because nothing in the design justifies a
                     directional test.

Every P value drawn or written out is computed here from the per-animal data; none is hard-coded.
Panel A is drawn only when --crop-dir is given, since the cropped PAS images are not redistributed.

Usage: plot_aging_phenotype.py --raw F [--crop-dir D] --out-prefix P
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tukey import tukey_hsd

ap = argparse.ArgumentParser()
ap.add_argument("--raw", required=True, help="source_data/Figure2_aging_cohort_raw.csv")
ap.add_argument("--crop-dir", default=None, help="directory with the cropped PAS images (panel A)")
ap.add_argument("--out-prefix", required=True)
a = ap.parse_args()

df = pd.read_csv(a.raw)
FULL = {"B6": "C57BL/6J", "FVB": "FVB/N"}
COL = {"B6": "#8172B2", "FVB": "#C44E52"}
STATS = []


def two_way_anova(sub, log=False):
    """strain x age two-way ANOVA; log=True fits the log-transformed response."""
    d = sub.copy()
    d["y"] = np.log(d.value) if log else d.value
    return sm.stats.anova_lm(smf.ols("y ~ C(strain) * C(age_label)", data=d).fit(), typ=2)


def holm(pvals):
    p = np.asarray(pvals, float); n = len(p); adj = np.empty(n); run = 0.0
    for rank, i in enumerate(np.argsort(p)):
        run = max(run, (n - rank) * p[i]); adj[i] = min(run, 1.0)
    return adj


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def trajectory(ax, measure, ages, ylabel, title, age_labels, test, anova_log=False, subtitle=None,
               title_pad=None):
    """Per-strain trajectory with individual animals and mean +/- SD.

    test: 'tukey'  Tukey-Kramer over all strain-by-age groups, FVB vs B6 read off within each age
          'holm'   Welch's t-test per age, Holm-corrected across the ages of this panel
          'welch'  Welch's t-test per age, uncorrected
    """
    sub = df[df.measure == measure]
    vals = {"B6": {}, "FVB": {}}
    rng = np.random.default_rng(0)
    for s in ["B6", "FVB"]:
        means, sds, xs = [], [], []
        for i, ag in enumerate(ages):
            v = sub[(sub.strain == s) & (sub.age_label == ag)].value
            vals[s][i] = v
            ax.scatter(rng.normal(i, 0.055, len(v)), v, s=22, color=COL[s],
                       alpha=.55, edgecolor="white", linewidth=.5, zorder=3)
            means.append(v.mean()); sds.append(v.std()); xs.append(i)
            STATS.append(dict(panel=title[0], measure=measure, group=f"{s} {ag}", quantity="mean",
                              value=v.mean()))
            STATS.append(dict(panel=title[0], measure=measure, group=f"{s} {ag}", quantity="SD",
                              value=v.std()))
        ax.errorbar(xs, means, yerr=sds, color=COL[s], lw=1.8, marker="o", ms=5,
                    capsize=3.5, zorder=4, label=FULL[s])

    if test == "tukey":
        g = sub.strain + "_" + sub.age_label
        tk = tukey_hsd(sub.value.values, g.values)
        p_used = [next(v for k, v in tk.items() if set(k) == {f"B6_{ag}", f"FVB_{ag}"}) for ag in ages]
        label = "Tukey-Kramer, six strain-by-age groups"
    else:
        raw = [stats.ttest_ind(vals["FVB"][i], vals["B6"][i], equal_var=False).pvalue
               for i in range(len(ages))]
        p_used = holm(raw) if test == "holm" else raw
        label = ("Welch t-test, Holm-corrected across the two ages" if test == "holm"
                 else "Welch t-test on untransformed values, uncorrected")
        for ag, pr in zip(ages, raw):
            STATS.append(dict(panel=title[0], measure=measure, group=f"FVB vs B6 {ag}",
                              quantity="Welch P (unadjusted)", value=pr))
    for ag, p in zip(ages, p_used):
        STATS.append(dict(panel=title[0], measure=measure, group=f"FVB vs B6 {ag}",
                          quantity=f"P plotted ({label})", value=p))

    ax.set_xticks(range(len(ages))); ax.set_xticklabels(age_labels, fontsize=11.2)
    ax.set_ylabel(ylabel, fontsize=12.6)
    ax.set_title(title, fontsize=13.4, weight="bold", loc="left", pad=title_pad)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11.2)
    ax.legend(fontsize=10.5, frameon=False)

    y0, y1 = ax.get_ylim(); pad = (y1 - y0) * 0.09
    for i, p in enumerate(p_used):
        top_i = max(vals["B6"][i].max(), vals["FVB"][i].max())
        ax.text(i, top_i + pad * 0.2, stars(p), ha="center", va="bottom",
                fontsize=11.9, weight="bold", clip_on=False)
    ax.set_ylim(y0, y1 + pad)
    if subtitle:
        ax.text(0, 1.03, subtitle(p_used), transform=ax.transAxes, ha="left", va="bottom",
                fontsize=9.2, color="#666666")

    aov = two_way_anova(sub, log=anova_log)
    p_s = aov.loc["C(strain)", "PR(>F)"]; p_a = aov.loc["C(age_label)", "PR(>F)"]
    p_i = aov.loc["C(strain):C(age_label)", "PR(>F)"]
    for term, key in [("strain", "C(strain)"), ("age", "C(age_label)"),
                      ("strain x age", "C(strain):C(age_label)")]:
        STATS.append(dict(panel=title[0], measure=measure, group=f"two-way ANOVA{' (log)' if anova_log else ''}",
                          quantity=f"{term} F({int(aov.loc[key,'df'])},{int(aov.loc['Residual','df'])})",
                          value=aov.loc[key, "F"]))
        STATS.append(dict(panel=title[0], measure=measure, group=f"two-way ANOVA{' (log)' if anova_log else ''}",
                          quantity=f"{term} P", value=aov.loc[key, "PR(>F)"]))
    ax.text(0, -0.30,
            f"two-way ANOVA{' (log)' if anova_log else ''}\n"
            f"strain $p$ = {p_s:.2g}\nage $p$ = {p_a:.2g}\nstrain × age $p$ = {p_i:.2g}",
            transform=ax.transAxes, ha="left", va="top", fontsize=9.5, color="#333333",
            linespacing=1.4,
            bbox=dict(boxstyle="round,pad=.3", fc="#f5f5f5", ec="#bbbbbb", lw=.5))


CROP_UM, BAR_UM = 150.0, 50.0
fig = plt.figure(figsize=(13.2, 9.3))
gs = fig.add_gridspec(2, 12, height_ratios=[1.0, 1.06], hspace=.26, wspace=2.2)

# --- A: representative PAS images (drawn only when the crops are available) ---
if a.crop_dir and os.path.isdir(a.crop_dir):
    panels = [("B6_PAS_1", "C57BL/6J", "10 wk"), ("B6_PAS_2", "C57BL/6J", "1 year"),
              ("FVB_PAS_1", "FVB/N", "10 wk"), ("FVB_PAS_2", "FVB/N", "1 year")]
    for i, (fname, strain, age) in enumerate(panels):
        ax = fig.add_subplot(gs[0, i * 3:(i + 1) * 3])
        ax.imshow(mpimg.imread(f"{a.crop_dir}/{fname}_crop.png"), extent=[0, CROP_UM, 0, CROP_UM])
        ax.set_xticks([]); ax.set_yticks([])
        c = COL["B6"] if strain == "C57BL/6J" else COL["FVB"]
        for sp in ax.spines.values():
            sp.set_edgecolor(c); sp.set_linewidth(1.6)
        ax.set_title(f"{strain}   {age}", fontsize=12.6, color=c, weight="bold", pad=4)
        x0, y0 = CROP_UM - BAR_UM - 8, 8
        ax.add_patch(Rectangle((x0, y0), BAR_UM, 2.6, color="black", zorder=5))
        if i == 0:
            ax.text(x0 + BAR_UM / 2, y0 + 4.5, f"{BAR_UM:.0f} µm", ha="center",
                    va="bottom", fontsize=10.5, zorder=5)
    fig.text(.013, .975, "A  Representative glomeruli (PAS)", fontsize=14.7,
             weight="bold", va="top", ha="left")

axB = fig.add_subplot(gs[1, 0:4]); axC = fig.add_subplot(gs[1, 4:8]); axD = fig.add_subplot(gs[1, 8:12])

trajectory(axB, "kidney_weight", ["8wk", "1yr", "1.5yr"], "Kidney weight (g)",
           "B  Kidney mass with age", ["8 wk", "1 year", "1.5 yr"], test="tukey")
trajectory(axC, "glom_diameter", ["10wk", "1yr"], "Glomerular long-axis diameter (µm)",
           "C  Glomerular long-axis\n    diameter", ["10 wk", "1 year"], test="holm",
           subtitle=lambda p: f"1 year: Holm-adjusted $P$ = {p[1]:.1e}".replace("e-0", " × 10$^{-") + "}$")
trajectory(axD, "uacr", ["10wk", "52wk"], "UACR (mg/g creatinine)",
           "D  Urinary albumin-to-\n    creatinine ratio", ["10 wk", "52 wk"], test="welch",
           anova_log=True)

# UACR sensitivity analyses: two-sided exact Mann-Whitney, and the within-strain age comparisons.
ua = df[df.measure == "uacr"]
for ag in ["10wk", "52wk"]:
    b6 = ua[(ua.strain == "B6") & (ua.age_label == ag)].value
    fv = ua[(ua.strain == "FVB") & (ua.age_label == ag)].value
    STATS.append(dict(panel="D", measure="uacr", group=f"FVB vs B6 {ag}",
                      quantity="Mann-Whitney U, two-sided exact P",
                      value=stats.mannwhitneyu(fv, b6, alternative="two-sided", method="exact").pvalue))
for s in ["FVB", "B6"]:
    v1 = ua[(ua.strain == s) & (ua.age_label == "10wk")].value
    v2 = ua[(ua.strain == s) & (ua.age_label == "52wk")].value
    STATS.append(dict(panel="D", measure="uacr", group=f"{FULL[s]} 10wk vs 52wk",
                      quantity="Welch P (untransformed, unadjusted)",
                      value=stats.ttest_ind(v2, v1, equal_var=False).pvalue))

fig.text(.055, .012,
         "Male mice only. B–D: points, individual animals; bars, mean ± SD. "
         "B: Tukey–Kramer (six groups). C: Welch t-test, Holm across the two ages. "
         "D: Welch t-test on untransformed values. "
         "A: all four images at identical magnification.",
         ha="left", va="bottom", fontsize=9.6, color="#555555")
fig.subplots_adjust(left=.055, right=.99, top=.89, bottom=.23)
for ext in ("png", "pdf"):
    fig.savefig(f"{a.out_prefix}.{ext}", dpi=300, bbox_inches="tight")
pd.DataFrame(STATS).to_csv(f"{a.out_prefix}_statistics.tsv", sep="\t", index=False)
print(f"saved {a.out_prefix}")
