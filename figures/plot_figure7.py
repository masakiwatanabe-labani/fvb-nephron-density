#!/usr/bin/env python3
"""New Figure 7: unrestricted gene-set analysis of the developing-kidney transcriptome.

v2 (2026-09-26): panel B now shows the two predefined MSigDB Hallmark interferon sets in full,
at both stages. The previous panel used 18 genes taken from the leading edge of the P1 result, which
made it post-hoc; those 18 genes are no longer shown in the main figure.

An earlier version of the analysis tested only gene sets whose NAME matched a kidney regular expression
(106 sets).
Removing that filter (10,082 sets) reveals two programmes at P1. Panel D reports their effect sizes
alongside canonical proliferation and growth-factor markers, so the weakness of the metabolic shift is
visible rather than implied.

Style follows the other manuscript figures (plot_figure7.py): no top/right spines, 11-13 pt labels.
Usage: plot_new_figure7.py --gsea-p1 F --gsea-e135 F --panelb F --panelc F --deseq-p1 F --id2name F --out-prefix P
"""
import argparse
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ap = argparse.ArgumentParser()
ap.add_argument("--gsea-p1", required=True); ap.add_argument("--gsea-e135", required=True)
ap.add_argument("--panelb", required=True, help="per-gene log2FC of the two Hallmark interferon sets")
ap.add_argument("--panelc", required=True)
ap.add_argument("--deseq-p1", required=True); ap.add_argument("--id2name", required=True)
ap.add_argument("--out-prefix", required=True)
a = ap.parse_args()

C_IFN, C_MET, C_NS, C_MARK = "#C44E52", "#4C72B0", "#C8C8C8", "#55A868"
def style(ax):
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=10.5)

IMMUNE = ("INTERFERON|VIRUS|VIRAL|ANTIGEN|MHC|IMMUN|DEFENSE|CYTOKINE|INNATE|TOLL|NFK|"
          "INFLAMMAT|LYMPHOCYTE|LEUKOCYTE|T_CELL|B_CELL|COMPLEMENT")

g = pd.read_csv(a.gsea_p1, sep="\t")
fig, axes = plt.subplots(1, 4, figsize=(21.8, 5.4),
                         gridspec_kw={"wspace": 0.34, "width_ratios": [1, 1, 1, 1.28]})

# ---------------- A: every tested set at P1 ----------------
ax = axes[0]; style(ax)
g = g.dropna(subset=["padj", "NES"]).copy()
g["y"] = -np.log10(g.padj.clip(lower=1e-16))
sig = g.padj < 0.05
imm = g.pathway.str.contains(IMMUNE)
ax.scatter(g.NES[~sig], g.y[~sig], s=7, c=C_NS, alpha=.55, linewidths=0, label=f"not significant ({int((~sig).sum())})")
ax.scatter(g.NES[sig & imm], g.y[sig & imm], s=17, c=C_IFN, linewidths=0, label=f"interferon / immune ({int((sig&imm).sum())})")
ax.scatter(g.NES[sig & ~imm], g.y[sig & ~imm], s=17, c=C_MET, linewidths=0, label=f"other ({int((sig&~imm).sum())})")
ax.axhline(-np.log10(0.05), color="k", lw=.8, ls=":")
ax.text(ax.get_xlim()[0], -np.log10(0.05), " FDR 0.05", va="bottom", ha="left", fontsize=9, color="k")
top = g.nsmallest(1, "padj").iloc[0]
ax.annotate("interferon gamma\nresponse", (top.NES, top.y), textcoords="offset points", xytext=(30, -34),
            fontsize=9, color=C_IFN, ha="left",
            arrowprops=dict(arrowstyle="-", color=C_IFN, lw=.7, shrinkA=2, shrinkB=3))
ax.set_xlabel("Normalised enrichment score (FVB/N vs C57BL/6J)", fontsize=11)
ax.set_ylabel("$-\\log_{10}$ adjusted $P$", fontsize=11)
ax.set_title(f"A  Without a kidney name filter,\n{len(g):,} sets tested at P1", fontsize=12.6, weight="bold")
ax.legend(fontsize=8.6, frameon=False, loc="center right", borderaxespad=.4)

# ---------------- B: predefined Hallmark interferon sets, E13.5 vs P1 ----------------
ax = axes[1]; style(ax)
b = pd.read_csv(a.panelb, sep="\t")
SETLAB = {"HALLMARK_INTERFERON_GAMMA_RESPONSE": "interferon\u2011\u03b3\nresponse",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE": "interferon\u2011\u03b1\nresponse"}
order = [(s_, tp) for s_ in SETLAB for tp in ("E13.5", "P1")]
data, labs, cols, notes = [], [], [], []
for s_, tp in order:
    sub = b[(b.set == s_) & (b.timepoint == tp)].dropna(subset=["log2FoldChange"])
    data.append(sub.log2FoldChange.values)
    labs.append(f"{SETLAB[s_]}\n{tp}")
    cols.append(C_NS if tp == "E13.5" else C_IFN)
    frac = (sub.padj < 0.1).sum() / sub.padj.notna().sum() if sub.padj.notna().any() else float("nan")
    notes.append((len(sub), float(np.median(sub.log2FoldChange)), frac))
bp = ax.boxplot(data, positions=range(len(data)), widths=.62, patch_artist=True,
                showfliers=False, medianprops=dict(color="k", lw=1.6))
for patch, c in zip(bp["boxes"], cols):
    patch.set_facecolor(c); patch.set_alpha(.55); patch.set_edgecolor("#555555"); patch.set_linewidth(.8)
rng_b = np.random.default_rng(1)
for i, v in enumerate(data):
    ax.scatter(i + rng_b.uniform(-.17, .17, len(v)), v, s=5, c="#444444", alpha=.35, linewidths=0, zorder=3)
ax.axhline(0, color="k", lw=.8)
ax.set_xticks(range(len(data))); ax.set_xticklabels(labs, fontsize=8.6)
ax.set_ylabel("log$_2$ fold change (FVB/N vs C57BL/6J)", fontsize=11)
ax.set_title("B  Predefined Hallmark interferon sets", fontsize=12.6, weight="bold", pad=20)
lo = min(np.percentile(v, 1) for v in data); hi = max(np.percentile(v, 99) for v in data)
pad = (hi - lo) * .34
ax.set_ylim(lo - pad * .5, hi + pad)
for i, (n_, med, frac) in enumerate(notes):
    ax.annotate(f"n = {n_}\nmed. {med:+.2f}",
                (i, hi + pad * .12), ha="center", va="bottom", fontsize=8.2, linespacing=1.25)
ax.text(.5, 1.015, "All available genes; no significance filter", transform=ax.transAxes,
        ha="center", va="bottom", fontsize=8.4, color="#777777")
ax.text(.5, -0.20, "Points: genes; boxes: median and interquartile range\n"
        "Whiskers: within 1.5 \u00d7 IQR; all values shown", transform=ax.transAxes,
        ha="center", va="top", fontsize=7.8, color="#777777")

# ---------------- C: kidney P1 vs adult liver ----------------
ax = axes[2]; style(ax)
c = pd.read_csv(a.panelc, sep="\t")
c["programme"] = c.programme.str.replace("\\n", "\n", regex=False)
# Omit programmes not supported by the final strand-specific quantification.
# Keep the retained programme medians exactly as supplied in the source table.
programme_key = c.programme.str.replace(r"\s+", " ", regex=True).str.strip().str.casefold()
c = c.loc[~programme_key.isin(["myc targets", "ribosome biogenesis"])].copy()
yy = np.arange(len(c))[::-1]
h = 0.36
ax.barh(yy + h/2, c.kidney_P1_median, height=h, color=[C_IFN if v < 0 else C_MET for v in c.kidney_P1_median],
        label="kidney, P1 (this study)")
ax.barh(yy - h/2, c.liver_median, height=h, color="white", edgecolor="#777777", linewidth=1.0,
        label="liver, adult (GSE167328)")
ax.axvline(0, color="k", lw=.8)
ax.set_yticks(yy); ax.set_yticklabels(c.programme, fontsize=8.8)
ax.set_xlabel("median log$_2$ fold change of leading-edge genes", fontsize=11)
ax.set_title("C  Interferon direction is reproduced\nin an independent liver dataset", fontsize=12.6, weight="bold")
ax.legend(fontsize=8.8, frameon=False, loc="lower right")

# ---------------- D: effect sizes, with proliferation markers ----------------
ax = axes[3]; style(ax)
gp = pd.read_csv(a.gsea_p1, sep="\t")
i2n = pd.read_csv(a.id2name, sep="\t", header=None, names=["gene_id", "name"])
d = pd.read_csv(a.deseq_p1, sep="\t").merge(i2n, on="gene_id")
d = d.sort_values(["name", "gene_id"]).drop_duplicates("name")
def le_of(names):
    s = set()
    for n in names:
        r = gp[gp.pathway == n]
        if len(r): s |= set(str(r.iloc[0].leadingEdge).split(";"))
    return s
ifn = d[d.name.isin(le_of(["HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE"]))]
met = d[d.name.isin(le_of(["HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM",
                           "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
                           "GOBP_MITOCHONDRION_ORGANIZATION"]))]
parts = ax.violinplot([ifn.log2FoldChange.dropna(), met.log2FoldChange.dropna()], positions=[0, 1.9],
                      widths=.78, showextrema=False, showmedians=True)
for pc, col in zip(parts["bodies"], [C_IFN, C_MET]):
    pc.set_facecolor(col); pc.set_alpha(.45); pc.set_edgecolor("none")
parts["cmedians"].set_color("k"); parts["cmedians"].set_linewidth(1.3)
MARK = ["Mki67", "Top2a", "Pcna", "Myc", "Igf1", "Mtor"]
mk = d.set_index("name").reindex(MARK).reset_index()
xm = 3.6 + np.arange(len(mk)) * 0.52
ax.scatter(xm, mk.log2FoldChange, s=44, c=C_MARK, edgecolors="k", linewidths=.5, zorder=3)
ax.axhline(0, color="k", lw=.8)
ax.set_xticks([0, 1.9] + list(xm))
ax.set_xticklabels([f"interferon\n\n(n = {len(ifn)})",
                    f"mitochondrial /\noxidative\n(n = {len(met)})"] + list(mk.name), fontsize=9.0)
for t in ax.get_xticklabels()[2:]:
    t.set_fontsize(8.6); t.set_style("italic"); t.set_rotation(45); t.set_ha("right")
ax.set_xlim(-.75, xm[-1] + .5)
ax.set_ylim(-5, 3)
ax.text(.99, .985, "tails beyond the axis not shown", transform=ax.transAxes,
        ha="right", va="top", fontsize=8.2, color="#666666")
ax.set_ylabel("log$_2$ fold change", fontsize=11)
ax.set_title("D  Metabolic shifts and selected\nproliferation-related transcripts", fontsize=12.6, weight="bold")
for xp, v in [(0, ifn.log2FoldChange.median()), (1.9, met.log2FoldChange.median())]:
    ax.annotate(f"{v:+.2f}", (xp, v), textcoords="offset points", xytext=(0, 7),
                ha="center", fontsize=8.8, weight="bold")

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{a.out_prefix}.{ext}", dpi=300, bbox_inches="tight")
vals = pd.DataFrame([
    {"panel": "A", "item": "sets tested at P1", "value": len(g)},
    {"panel": "A", "item": "sets padj<0.05", "value": int(sig.sum())},
    {"panel": "A", "item": "of which interferon/immune", "value": int((sig & imm).sum())},
    {"panel": "A", "item": "minimum padj", "value": float(g.padj.min())},
] + [{"panel": "B", "item": f"{s_.replace('HALLMARK_','')} {tp} {k}", "value": v}
      for (s_, tp), (n_, med, frac) in zip(order, notes)
      for k, v in [("n", n_), ("median log2FC", round(med, 4)), ("frac padj<0.1", round(frac, 3))]] + [
    {"panel": "D", "item": "interferon leading edge median log2FC", "value": float(ifn.log2FoldChange.median())},
    {"panel": "D", "item": "mito/oxidative leading edge median log2FC", "value": float(met.log2FoldChange.median())},
] + [{"panel": "D", "item": f"{r['name']} log2FC / padj",
      "value": f"{r.log2FoldChange:+.2f} / {r.padj:.3g}"} for _, r in mk.iterrows()])
vals.to_csv(f"{a.out_prefix}_plotted_values.tsv", sep="\t", index=False)
print(f"saved {a.out_prefix}")
