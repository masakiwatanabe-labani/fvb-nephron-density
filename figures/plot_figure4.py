#!/usr/bin/env python3
"""Figure 4 (differential expression and GSEA). Same layout and styling as the 'Figure 3' block of
plot_omics.py (drawn as Figure4 by an earlier version of the analysis), with paths as arguments
(no /tmp) and panel D's title and FDR note derived from the data.

v2 (2026-09-28): panels A and B additionally report how many genes were tested, how many received a
non-missing adjusted P value after DESeq2's independent filtering, and how many reached FDR < 0.1
before the fold-change cut. Those three counts and the annotated count are different quantities and
had been reported inconsistently; all four are now printed from the same DESeq2 table."""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ap = argparse.ArgumentParser()
ap.add_argument("--deseq-dir", required=True); ap.add_argument("--fgsea-tsv", required=True, help="fgsea_E13.5.tsv")
ap.add_argument("--id2name", required=True); ap.add_argument("--gsea-note", default=""); ap.add_argument("--out-prefix", required=True)
a = ap.parse_args()
C_FVB = "#C44E52"
def style(ax): ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=11.2)
fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.8))
vals_out = []
for j, tp in enumerate(["E13.5", "P1"]):
    raw = pd.read_csv(os.path.join(a.deseq_dir, f"DESeq2_{tp}.tsv"), sep="\t")
    n_tested = len(raw)                                   # genes that entered DESeq2
    d = raw.dropna(subset=["padj", "log2FoldChange"])
    n_padj = len(d)                                       # survived independent filtering
    n_fdr = int((d.padj < 0.1).sum())                     # FDR < 0.1, before the fold-change cut
    ax = axes[0, j]
    sig = (d.padj < 0.1) & (d.log2FoldChange.abs() > 0.5)
    ax.scatter(d.log2FoldChange[~sig], -np.log10(d.padj[~sig]), s=2.5, color="#cccccc", alpha=.5, rasterized=True)
    ax.scatter(d.log2FoldChange[sig], -np.log10(d.padj[sig]), s=4, color=C_FVB, alpha=.55, rasterized=True)
    ax.axhline(1, ls="--", lw=.7, color="#888"); ax.axvline(0, lw=.6, color="#888"); ax.set_xlim(-4, 4)
    ax.set_xlabel("log2FC (FVB vs B6)", fontsize=12.6); ax.set_ylabel("−log10 FDR", fontsize=12.6)
    ax.set_title(f"{'A' if j == 0 else 'B'}  Differential expression: {tp}", fontsize=14, weight="bold"); style(ax)
    ax.text(.02, .96, f"{int(sig.sum())} genes\nFDR<0.1, |log2FC|>0.5", transform=ax.transAxes,
            va="top", fontsize=9.8, color="#444")
    ax.text(.02, .845, f"Tested: {n_tested:,}; padj available: {n_padj:,}\nFDR < 0.1 only: {n_fdr:,}",
            transform=ax.transAxes, va="top", fontsize=8.6, color="#777")
    L = "A" if j == 0 else "B"
    vals_out += [(L, f"{tp} genes FDR<0.1 & |log2FC|>0.5", int(sig.sum())),
                 (L, f"{tp} genes tested", n_tested),
                 (L, f"{tp} genes with non-missing padj", n_padj),
                 (L, f"{tp} genes FDR<0.1 (no fold-change cut)", n_fdr)]
ax = axes[1, 0]
genes = ["Gdnf", "Gfra1", "Ret", "Six2", "Cited1", "Wnt9b", "Wnt11", "Etv4"]
id2n = pd.read_csv(a.id2name, sep="\t", header=None, names=["gid", "name"]); n2i = dict(zip(id2n.name, id2n.gid))
w = 0.36
for i, tp in enumerate(["E13.5", "P1"]):
    d = pd.read_csv(os.path.join(a.deseq_dir, f"DESeq2_{tp}.tsv"), sep="\t").set_index("gene_id")
    vals, sigs = [], []
    for g in genes:
        gid = n2i.get(g)
        if gid in d.index: vals.append(d.loc[gid, "log2FoldChange"]); sigs.append(d.loc[gid, "padj"] < 0.1)
        else: vals.append(np.nan); sigs.append(False)
        vals_out.append(("C", f"{g} {tp} log2FC / FDR<0.1", f"{vals[-1]:.3f} / {sigs[-1]}"))
    xs = np.arange(len(genes)) + (i - 0.5) * w
    ax.bar(xs, vals, w, color=["#4C72B0", "#55A868"][i], alpha=.85, label=tp)
    for x, v, s in zip(xs, vals, sigs):
        if s and not np.isnan(v): ax.text(x, v + (0.04 if v > 0 else -0.09), "*", ha="center", fontsize=15.4, weight="bold")
ax.axhline(0, color="black", lw=.8)
ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, fontsize=11.2, style="italic", rotation=30)
ax.set_ylabel("log2FC (FVB vs B6)", fontsize=12.6)
ax.set_title("C  GDNF/RET axis is not coordinately reduced", fontsize=14, weight="bold")
ax.legend(fontsize=11.2, frameon=False, loc="upper right", bbox_to_anchor=(1.06, 1.0)); style(ax)
ax.text(.98, .04, "* FDR<0.1", transform=ax.transAxes, ha="right", fontsize=9.8, color="#555")
ax = axes[1, 1]
allf = pd.read_csv(a.fgsea_tsv, sep="\t")
f = allf.sort_values("pval").head(8).copy()
f["short"] = f.pathway.str.replace("GOBP_", "").str.replace("REACTOME_", "").str.replace("_", " ").str.title().str[:34]
ax.barh(range(len(f)), f.NES, color=[C_FVB if n < 0 else "#4C72B0" for n in f.NES], alpha=.8)
ax.set_yticks(range(len(f))); ax.set_yticklabels(f.short, fontsize=9.8); ax.set_ylim(len(f) + 0.35, -0.6); ax.axvline(0, color="black", lw=.8)
ax.set_xlabel("NES", fontsize=12.6)
nsig = int((allf.padj < 0.1).sum())
title = "D  GSEA (E13.5): no set reaches\nFDR<0.1" if nsig == 0 else f"D  GSEA (E13.5): {nsig} set{'s' if nsig > 1 else ''} at\nFDR<0.1"
ax.set_title(title, fontsize=14, weight="bold"); style(ax)
ax.text(.02, .025, f"min FDR = {allf.padj.min():.2f} ({len(allf)} sets tested{a.gsea_note})", transform=ax.transAxes, fontsize=9.8, color="#555")
vals_out += [("D", "sets tested", len(allf)), ("D", "sets FDR<0.1", nsig), ("D", "min FDR", round(float(allf.padj.min()), 4))] + [("D", f"NES {p}", round(n, 3)) for p, n in zip(f.pathway, f.NES)]
fig.tight_layout(); fig.savefig(a.out_prefix + ".png", dpi=600, bbox_inches="tight"); fig.savefig(a.out_prefix + ".pdf", bbox_inches="tight"); plt.close(fig)
pd.DataFrame(vals_out, columns=["panel", "item", "value"]).to_csv(a.out_prefix + "_plotted_values.tsv", sep="\t", index=False)
print("saved", a.out_prefix)
