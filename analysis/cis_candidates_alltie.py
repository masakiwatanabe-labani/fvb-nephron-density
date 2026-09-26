#!/usr/bin/env python3
"""Candidate prioritisation keeping EVERY tied nearest-gene assignment.

cis_candidates_param.py:20 drops duplicates on (chr, start, var_id), which keeps only the first of
the genes that bedtools closest returned at equal distance. 295 of the 2,216 regulatory variants lie
inside two or more overlapping gene annotations (all at distance 0), and 316 gene assignments are
discarded by that rule, the survivor depending on input order. Here every tied (variant, gene) pair
is kept, so a variant inside overlapping genes is credited to all of them.

Peak membership is also tested per variant RECORD (as in rerun4/exact_match/cis_candidates_exact.py)
rather than per position, so CHROM:POS:REF:ALT identifies a variant throughout.

Everything else - the +-250 bp window, the DESeq2 thresholds and the QC - is unchanged.
Usage: same arguments as cis_candidates_param.py, plus --min-basemean.
"""
import argparse, os, subprocess
import numpy as np, pandas as pd
ap = argparse.ArgumentParser()
ap.add_argument("--reg-snps", required=True); ap.add_argument("--reg-indels", required=True)
ap.add_argument("--chip-dir", required=True); ap.add_argument("--deseq-dir", required=True)
ap.add_argument("--id2name", required=True); ap.add_argument("--matdir", required=True)
ap.add_argument("--gene-chr", default=None); ap.add_argument("--exclude-chroms", default="")
ap.add_argument("--min-basemean", type=float, default=25.0); ap.add_argument("--outdir", required=True)
a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True); W = os.path.join(a.outdir, "work"); os.makedirs(W, exist_ok=True)

cols = ["chr","start","end","var_id","g_chr","g_start","g_end","gene_id","dist"]
reg = pd.concat([pd.read_csv(a.reg_snps, sep="\t", header=None, names=cols),
                 pd.read_csv(a.reg_indels, sep="\t", header=None, names=cols)], ignore_index=True)
reg = reg.drop_duplicates(subset=["chr","start","var_id","gene_id"]); reg = reg[reg["dist"] <= 100000].copy()
CHIP = {"Six2":"Six2.grcm39.summit.bed","Six2_BF":"Six2_BF.grcm39.summit.bed","Osr1":"Osr1_BF.grcm39.summit.bed",
        "Wt1":"Wt1.grcm39.summit.bed","Hoxd11":"Hoxd11_BF.grcm39.summit.bed"}
# keep var_id in the BED so each record is tested on its own interval
reg[["chr","start","end","var_id"]].drop_duplicates().sort_values(["chr","start"]).to_csv(
    f"{W}/cis_regvars.bed", sep="\t", header=False, index=False)
hits = {}
for tf, fn in CHIP.items():
    cmd = (f"awk 'BEGIN{{OFS=\"\\t\"}} {{s=$2-250; if(s<0)s=0; print $1,s,$3+250}}' {os.path.join(a.chip_dir, fn)} "
           f"| sort -k1,1 -k2,2n | bedtools merge -i - > {W}/cis_peaks_{tf}.bed && "
           f"bedtools intersect -a {W}/cis_regvars.bed -b {W}/cis_peaks_{tf}.bed -u > {W}/cis_in_{tf}.bed")
    subprocess.run(["bash","-c",cmd], check=True)
    h = pd.read_csv(f"{W}/cis_in_{tf}.bed", sep="\t", header=None, names=["chr","start","end","var_id"])
    hits[tf] = set(h["var_id"])                      # <- keyed by record, not by position
reg["chip_tfs"] = reg["var_id"].map(lambda v: ";".join(t for t, s in hits.items() if v in s))
reg["in_chip"] = reg["chip_tfs"] != ""
agg = reg.groupby("gene_id").agg(n_reg_vars=("var_id","count"), n_in_chip=("in_chip","sum"),
    chip_tfs=("chip_tfs", lambda x: ";".join(sorted({t for v in x for t in v.split(";") if t}))),
    min_dist=("dist","min")).reset_index()
id2 = pd.read_csv(a.id2name, sep="\t", header=None, names=["gene_id","gene_name"])
agg["gene_name"] = agg["gene_id"].map(dict(zip(id2.gene_id, id2.gene_name)))
for tag, tp in [("E135","E13.5"), ("P1","P1")]:
    d = pd.read_csv(f"{a.deseq_dir}/DESeq2_{tp}.tsv", sep="\t").set_index("gene_id")
    for c, src in [("log2FC","log2FoldChange"), ("padj","padj"), ("baseMean","baseMean")]:
        agg[f"{c}_{tag}"] = d[src].reindex(agg["gene_id"]).values
def best(row):
    c = [(row[f"padj_{t}"], abs(row[f"log2FC_{t}"]), t) for t in ["E135","P1"]
         if pd.notna(row[f"padj_{t}"]) and pd.notna(row[f"log2FC_{t}"])]
    if not c: return pd.Series({"best_tp": None, "best_padj": None, "best_absLFC": None})
    p, l, t = min(c); return pd.Series({"best_tp": t, "best_padj": p, "best_absLFC": l})
agg = pd.concat([agg, agg.apply(best, axis=1)], axis=1); agg.to_csv(f"{a.outdir}/cis_discovery_all.tsv", sep="\t", index=False)
sig = agg[(agg["best_padj"] < 0.1) & (agg["best_absLFC"] > 0.5)]
cand = sig[sig["n_in_chip"] > 0].sort_values("best_padj"); cand.to_csv(f"{a.outdir}/cis_discovery_candidates.tsv", sep="\t", index=False)
excl = [c for c in a.exclude_chroms.split(",") if c]
gchr = pd.read_csv(a.gene_chr, sep="\t", usecols=["gene_id","chr"]).set_index("gene_id")["chr"] if a.gene_chr else None
def loadm(n):
    m = pd.read_csv(f"{a.matdir}/{n}_counts.tsv", sep="\t", index_col=0)
    return m[~gchr.reindex(m.index).isin(excl).values] if excl else m
mats = {n: loadm(n) for n in ["uncorrected","condA","condB"]}
def grp(mat, tp):
    suf = "135" if tp == "E135" else "P1"
    return [c for c in mat.columns if c.startswith("B6") and suf in c], [c for c in mat.columns if c.startswith("FVB") and suf in c]
def cpm_log2fc(mat, gid, tp):
    if gid not in mat.index: return np.nan
    b6, fvb = grp(mat, tp); cpm = mat.loc[gid] / mat.sum(axis=0) * 1e6
    mb, mf = cpm[b6].mean(), cpm[fvb].mean()
    return np.nan if mb <= 0 or mf <= 0 else np.log2(mf / mb)
def separation(mat, gid, tp):
    if gid not in mat.index: return False
    b6, fvb = grp(mat, tp); cpm = mat.loc[gid] / mat.sum(axis=0) * 1e6
    return bool(cpm[fvb].max() < cpm[b6].min() or cpm[fvb].min() > cpm[b6].max())
rows = []
for _, r in cand.iterrows():
    gid, tp = r["gene_id"], r["best_tp"]; lfc = {n: cpm_log2fc(m, gid, tp) for n, m in mats.items()}
    shift = abs(lfc["uncorrected"] - lfc["condB"]) if not (np.isnan(lfc["uncorrected"]) or np.isnan(lfc["condB"])) else np.nan
    rows.append({"gene_name": r["gene_name"], "gene_id": gid, "best_tp": tp, "baseMean": r[f"baseMean_{tp}"],
                 "log2FC_deseq": r[f"log2FC_{tp}"], "padj": r[f"padj_{tp}"], "lfc_uncorr": lfc["uncorrected"],
                 "lfc_condA": lfc["condA"], "lfc_condB": lfc["condB"], "bias_shift": shift,
                 "separated": separation(mats["condB"], gid, tp), "n_reg_vars": r["n_reg_vars"],
                 "n_in_chip": r["n_in_chip"], "chip_tfs": r["chip_tfs"], "min_dist": r["min_dist"]})
qc = pd.DataFrame(rows); qc.to_csv(f"{a.outdir}/cis_discovery_qc.tsv", sep="\t", index=False)
passed = qc[(qc["baseMean"] >= a.min_basemean) & (qc["bias_shift"] < 0.2) & (qc["separated"]) & (qc["log2FC_deseq"].abs() < 5)]
passed.to_csv(f"{a.outdir}/cis_discovery_passed.tsv", sep="\t", index=False)
print(f"regulatory-variant genes={len(agg)}; DE={len(sig)}; with ChIP={len(cand)}; PASSED={len(passed)}")
print("candidates:", sorted(passed.gene_name))
