#!/usr/bin/env python3
"""Table 2 rebuilt from the final re-quantification (-p --countReadPairs -s 2, all three arms).

The submitted Table 2 mixed two vintages: mapping rates from the original STAR runs of every arm, and
"Assigned to genes" from featureCounts run with -p only (read-level, strand-unaware). Here the mapping
rates come from the STAR runs that produced the counts actually used, and Assigned is the
fragment-level rate that featureCounts reports for the same run.

Usage: build_table2_final.py --samplesheet F --star-uncorrected DIR --star-rerun DIR
                             --counts-uncorrected DIR --counts-rerun DIR --fastp DIR --out F
"""
import argparse, os, re, json, glob
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--samplesheet", required=True)
ap.add_argument("--star-uncorrected", required=True, help="results/phase2/star_uncorrected")
ap.add_argument("--star-rerun", required=True, help="rerun3/S2/out/phase2_s2pairs/star (deleted after counting)")
ap.add_argument("--counts-uncorrected", required=True); ap.add_argument("--counts-rerun", required=True)
ap.add_argument("--star-fallback", default=None, help="prefix such as results/phase2/star_cond");
ap.add_argument("--fastp", default=None); ap.add_argument("--out", required=True)
a = ap.parse_args()

def star_log(path):
    if not os.path.exists(path): return {}
    d = {}
    for line in open(path):
        if "|" in line:
            k, v = line.split("|", 1); d[k.strip()] = v.strip()
    return d

def fc_summary(path):
    if not os.path.exists(path): return None
    s = pd.read_csv(path, sep="\t", index_col=0).iloc[:, 0]
    tot = s.sum()
    return 100 * s.get("Assigned", 0) / tot if tot else None

meta = pd.read_csv(a.samplesheet, sep="\t")
col = "sample" if "sample" in meta.columns else meta.columns[0]
rows = []
for _, m in meta.iterrows():
    s = str(m[col])
    u = star_log(f"{a.star_uncorrected}/{s}/Log.final.out")
    rec = {"Sample": s,
           "Strain": m.get("strain", ""), "Stage": m.get("stage", ""), "Sex": m.get("sex", ""),
           "Input read pairs (M)": round(int(u.get("Number of input reads", 0)) / 1e6, 2) if u else None,
           "Uniquely mapped, reference (%)": float(u["Uniquely mapped reads %"].rstrip("%")) if u else None}
    for arm, lab in [("condA", "cond. A"), ("condB", "cond. B")]:
        d = star_log(f"{a.star_rerun}/{s}.{arm}/Log.final.out")
        src = "rerun"
        if not d and a.star_fallback:
            # the rerun deleted each STAR directory together with its BAM, so its Log.final.out is gone.
            # The earlier run used the same STAR build, index, trimmed FASTQ and parameters
            # (run_phase2_continuation.sh:18-27 vs run_phase2_s2pairs.sh:30-37); its log is used here and
            # the provenance is recorded per row.
            d = star_log(f"{a.star_fallback}{arm[-1]}/{s}/Log.final.out"); src = "earlier run (identical parameters)"
        rec[f"Uniquely mapped, {lab} (%)"] = float(d["Uniquely mapped reads %"].rstrip("%")) if d else None
        rec[f"_source, {lab}"] = src if d else "NOT AVAILABLE"
    rec["Assigned fragments, reference (%)"] = fc_summary(f"{a.counts_uncorrected}/{s}.uncorrected.s2p.counts.txt.summary")
    for arm, lab in [("condA", "cond. A"), ("condB", "cond. B")]:
        rec[f"Assigned fragments, {lab} (%)"] = fc_summary(f"{a.counts_rerun}/{s}.{arm}.counts.txt.summary")
    rows.append(rec)
t = pd.DataFrame(rows)
for c in t.columns:
    if t[c].dtype.kind == "f": t[c] = t[c].round(2)
t.to_csv(a.out, index=False)
print(t.to_string(index=False))
missing = [c for c in t.columns if t[c].isna().all()]
if missing: print(f"\nCOLUMNS WITH NO DATA (source files absent): {missing}")
