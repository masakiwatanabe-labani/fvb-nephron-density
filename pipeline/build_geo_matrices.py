#!/usr/bin/env python3
"""Assemble GEO processed-data matrices from featureCounts per-sample outputs.

One file per processing arm, each a gene x sample table with annotation columns
(gene_id, gene_name, chr, length) followed by the twelve samples. TPM is computed from the same
counts and the featureCounts effective gene length (union of exons), so counts and TPM are consistent.

Usage: build_geo_matrices.py --counts-dir DIR --pattern 'SAMPLE.<tag>.counts.txt' --id2name F
                             --out-prefix P [--samples a,b,...]
"""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--counts-dir", required=True)
ap.add_argument("--pattern", required=True, help="file name with SAMPLE where the sample name goes")
ap.add_argument("--id2name", required=True)
ap.add_argument("--out-prefix", required=True)
ap.add_argument("--samples", default="B6135-1,B6135-2,B6135-3,B6P1-1,B6P1-2,B6P1-3,"
                                     "FVB135-1,FVB135-2,FVB135-3,FVBP1-1,FVBP1-2,FVBP1-3")
a = ap.parse_args()
samples = a.samples.split(",")

ann, cols = None, {}
for s in samples:
    f = os.path.join(a.counts_dir, a.pattern.replace("SAMPLE", s))
    d = pd.read_csv(f, sep="\t", comment="#")
    if ann is None:
        ann = d[["Geneid", "Chr", "Length"]].copy()
        # featureCounts lists one entry per exon; every exon of a gene is on the same contig here
        ann["Chr"] = ann["Chr"].astype(str).str.split(";").str[0]
        ann.columns = ["gene_id", "chr", "length"]
    else:
        assert (d["Geneid"].values == ann["gene_id"].values).all(), f"gene order differs in {f}"
    cols[s] = d.iloc[:, -1].values

m = pd.DataFrame(cols, index=ann["gene_id"].values)
name = pd.read_csv(a.id2name, sep="\t", header=None, names=["gene_id", "gene_name"]).drop_duplicates("gene_id")
ann = ann.merge(name, on="gene_id", how="left")
ann["gene_name"] = ann["gene_name"].fillna("")
out = pd.concat([ann.set_index("gene_id")[["gene_name", "chr", "length"]], m], axis=1)
out.index.name = "gene_id"
out.to_csv(a.out_prefix + "_counts.tsv", sep="\t")

rate = m.div(ann.set_index("gene_id")["length"].reindex(m.index), axis=0)
tpm = (rate / rate.sum() * 1e6).round(4)
outt = pd.concat([ann.set_index("gene_id")[["gene_name", "chr", "length"]], tpm], axis=1)
outt.index.name = "gene_id"
outt.to_csv(a.out_prefix + "_TPM.tsv", sep="\t")
print(f"{a.out_prefix}: {out.shape[0]} genes x {len(samples)} samples; "
      f"total counts (M) {(m.sum()/1e6).round(1).to_dict()}")
