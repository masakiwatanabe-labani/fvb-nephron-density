#!/usr/bin/env python3
"""Table 3 rebuild with a corrected summit-distance column, keeping every tied nearest-gene pair.

v3 of the builder. In addition to the two distance defects fixed in v2, the de-duplication on
(chr, start, var_id) is replaced by (chr, start, var_id, gene_id): the former kept only the first of
the genes that bedtools closest returned at equal distance, so a gene whose variants were all
assigned to an overlapping neighbour had no rows here and could not appear in the table.

Two defects in build_table3_param.py are fixed here:
  1. the summit dictionary omitted Six2_BF, so a variant qualified by a Six2(BioTagFLAG) peak had its
     distance measured against the conventional Six2 summit set - a different dataset. This is why
     Trnp1 (265 bp) and Usp29 (259 bp) appeared to breach the +-250 bp selection window.
  2. candidate summits were searched only within +-400 bp of the variant, which can return a summit
     farther away than the one that qualified the variant, and yields a blank when none is in range.
The distance is now the minimum over exactly the TF sets whose +-250 bp window contains the variant,
so the reported distance always satisfies the stated criterion and names the dataset that qualified it.
"""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--candidates-dir", required=True, help="dir with cis_discovery_passed.tsv / cis_discovery_candidates.tsv")
ap.add_argument("--reg-snps", required=True); ap.add_argument("--reg-indels", required=True)
ap.add_argument("--cis-in-dir", required=True, help="dir with cis_in_<TF>.bed")
ap.add_argument("--summit-dir", required=True)
ap.add_argument("--window", type=int, default=250, help="selection half-window around a summit (bp)")
ap.add_argument("--out", required=True)
a = ap.parse_args()

# summit file per TF key, including Six2_BF (the omission fixed here)
SUMMIT = {"Six2": "Six2", "Six2_BF": "Six2_BF", "Osr1": "Osr1_BF", "Wt1": "Wt1", "Hoxd11": "Hoxd11_BF"}
LABEL = {"Six2_BF": "Six2(BioTagFLAG)"}

p = pd.read_csv(f"{a.candidates_dir}/cis_discovery_passed.tsv", sep="\t")
c = pd.read_csv(f"{a.candidates_dir}/cis_discovery_candidates.tsv", sep="\t").set_index("gene_id")
cols = ["chr","start","end","var_id","g_chr","g_start","g_end","gene_id","dist"]
reg = pd.concat([pd.read_csv(a.reg_snps, sep="\t", header=None, names=cols),
                 pd.read_csv(a.reg_indels, sep="\t", header=None, names=cols)],
                ignore_index=True).drop_duplicates(subset=["chr","start","var_id","gene_id"])

chip = {}
for tf in SUMMIT:
    h = pd.read_csv(f"{a.cis_in_dir}/cis_in_{tf}.bed", sep="\t", header=None, names=["chr","start","end"])
    chip[tf] = set(zip(h.chr, h.start))
reg["tfs"] = reg.apply(lambda r: ";".join(tf for tf, s in chip.items() if (r.chr, r.start) in s), axis=1)
reg_chip = reg[reg.tfs != ""]

summits = {tf: pd.read_csv(f"{a.summit_dir}/{f}.grcm39.summit.bed", sep="\t", header=None,
                           names=["c","s","e","id","score"]) for tf, f in SUMMIT.items()}

def nearest(chrom, pos, tfs):
    """Minimum distance to a summit, over exactly the TF sets that qualified this variant."""
    best = (None, None)
    for tf in tfs:
        sub = summits[tf]
        sub = sub[sub.c == chrom]
        if len(sub) == 0: continue
        d = int((sub.e - pos).abs().min())
        if best[0] is None or d < best[0]: best = (d, tf)
    return best

rows, audit = [], []
for _, r in p.iterrows():
    v = reg_chip[reg_chip.gene_id == r.gene_id]
    if len(v) == 0: continue
    v0 = v.iloc[0]
    qual = [t for t in v0.tfs.split(";") if t]
    dmin, tfmin = nearest(v0.chr, v0.end, qual)
    # the old column, reproduced, for the old->new table
    old_sets = ["Six2","Osr1","Wt1","Hoxd11"]
    od, ot = None, None
    for tf in old_sets:
        sub = summits[tf]
        sub = sub[(sub.c == v0.chr) & (sub.e.between(v0.end - 400, v0.end + 400))]
        for _, x in sub.iterrows():
            dd = abs(int(x.e) - int(v0.end))
            if od is None or dd < od: od, ot = dd, tf
    audit.append({"Gene": r.gene_name, "Variant": f"{v0.chr}:{v0.end}", "qualifying_TFs": v0.tfs,
                  "old_distance": od, "old_factor": ot, "new_distance": dmin, "new_factor": tfmin,
                  "within_window": (dmin is not None and dmin <= a.window)})
    cc = c.loc[r.gene_id]
    rows.append({
        "Gene": r.gene_name, "Ensembl ID": r.gene_id, "Variant": f"{v0.chr}:{v0.end}",
        "TF peaks containing variant": ";".join(LABEL.get(t, t) for t in qual),
        "Distance to nearest summit (bp)": dmin,
        "Nearest summit factor": LABEL.get(tfmin, tfmin),
        "log2FC E13.5": round(cc.log2FC_E135, 3) if pd.notna(cc.log2FC_E135) else None,
        "FDR E13.5": f"{cc.padj_E135:.2e}" if pd.notna(cc.padj_E135) else None,
        "log2FC P1": round(cc.log2FC_P1, 3) if pd.notna(cc.log2FC_P1) else None,
        "FDR P1": f"{cc.padj_P1:.2e}" if pd.notna(cc.padj_P1) else None,
        "Mean normalised count": round(r.baseMean, 0),
    })
t3 = pd.DataFrame(rows)
t3["_s"] = t3[["FDR E13.5","FDR P1"]].apply(lambda x: min([float(v) for v in x if v is not None] or [1]), axis=1)
t3 = t3.sort_values("_s").drop(columns="_s")
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
t3.to_csv(a.out, index=False)
ad = pd.DataFrame(audit)
ad.to_csv(a.out.replace(".csv", "_distance_audit.tsv"), sep="\t", index=False)
breach = ad[~ad.within_window]
print(f"Table 3: {len(t3)} genes; distances recomputed over qualifying TF sets only")
print(f"rows whose reported distance changed: {int((ad.old_distance != ad.new_distance).sum())}")
print(f"rows still outside +-{a.window} bp: {len(breach)}")
if len(breach): print(breach.to_string(index=False))
