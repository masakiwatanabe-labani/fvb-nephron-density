#!/usr/bin/env python3
"""allreg_full.tsv: regulatory SNVs with columns c, s, e, id, pos, six2_pred.
six2_pred = the SNV position (seqnames, start of the motifbreakR GRanges) has >=1 motifbreakR disruption call for SIX2."""
import argparse, pandas as pd
ap = argparse.ArgumentParser(); ap.add_argument("--reg", required=True); ap.add_argument("--motifbreakr", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
r = pd.read_csv(a.reg, sep="\t", header=None, names=["c", "s", "e", "id"])
mb = pd.read_csv(a.motifbreakr, sep="\t", usecols=["seqnames", "start", "geneSymbol"])
six2 = set(zip(mb[mb.geneSymbol.str.upper() == "SIX2"].seqnames, mb[mb.geneSymbol.str.upper() == "SIX2"].start))
r["pos"] = r.e; r["six2_pred"] = [(c, p) in six2 for c, p in zip(r.c, r.pos)]
r.to_csv(a.out, sep="\t", index=False)
