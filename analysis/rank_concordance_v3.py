#!/usr/bin/env python3
"""Four-strain allele-dosage concordance on exactly the records of the allele-wise private set.

Denominator = records of the set defined in Results 4.6 that have a complete genotype in all four
phenotyped strains. Records are matched on CHROM/POS/REF/ALT, so records that merely overlap a
private indel are excluded (the previous run matched on position ranges and admitted 3,118 such
positions).

Strain ranking: glomerular density per gram of kidney weight, C57BL/6J > BALB/cJ > DBA/2J > FVB/N.
Concordant = dosage non-decreasing as density decreases, with at least one real increase. C57BL/6J
dosage is fixed at 0, which is what the reference genome implies.

Usage: rank_concordance_v3.py --geno-dir DIR --priv-snps VCF --priv-indels VCF --candidates TSV
                              --reg-snps BED --reg-indels BED --outdir DIR
"""
import argparse, os, subprocess
import pandas as pd
from scipy.stats import binomtest

ap = argparse.ArgumentParser()
for k in ["geno-dir", "priv-snps", "priv-indels", "candidates", "reg-snps", "reg-indels", "outdir"]:
    ap.add_argument(f"--{k}", required=True)
a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)

def dosage(gt):
    gt = gt.strip()
    if gt in (".", "./.", ".|."): return None
    return sum(1 for x in gt.replace("|", "/").split("/") if x not in ("0", "."))

def priv_records(vcf):
    out = subprocess.run(["bash", "-c", f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf}"],
                         capture_output=True, text=True).stdout
    return set(tuple(l.split("\t")) for l in out.strip().split("\n"))

priv = {"snps": priv_records(a.priv_snps), "indels": priv_records(a.priv_indels)}
frames, audit = [], []
for v in ["snps", "indels"]:
    g = pd.read_csv(f"{a.geno_dir}/{v}_genotypes_v3.txt", sep="\t", header=None,
                    names=["chrom", "pos", "ref", "alt", "gts"], dtype={"chrom": str, "pos": str})
    n_rows = len(g)
    g["key"] = list(zip(g.chrom, g.pos, g.ref, g.alt))
    g["in_priv"] = [k in priv[v] for k in g.key]
    d = g.gts.str.rstrip(",").str.split(",")
    g["dos"] = d.apply(lambda x: [dosage(y) for y in x] if len(x) == 4 else None)
    g["complete"] = g.dos.apply(lambda x: x is not None and all(y is not None for y in x))
    audit.append(dict(kind=v, private_records=len(priv[v]), rows_returned=n_rows,
                      rows_matching_private=int(g.in_priv.sum()),
                      rows_not_private=int((~g.in_priv).sum()),
                      private_with_complete_gt=int((g.in_priv & g.complete).sum()),
                      private_without_complete_gt=int((g.in_priv & ~g.complete).sum())))
    frames.append(g[g.in_priv & g.complete])

allv = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
allv[["d_BALB", "d_DBA", "d_B6N", "d_FVB"]] = pd.DataFrame(allv.dos.tolist(), index=allv.index)
# kidney-mass ranking, C57BL/6J fixed at 0
allv["concordant"] = (0 <= allv.d_BALB) & (allv.d_BALB <= allv.d_DBA) & (allv.d_DBA <= allv.d_FVB) & (allv.d_FVB > 0)
allv.drop(columns=["dos", "key", "gts", "in_priv", "complete"]).to_csv(
    f"{a.outdir}/background_records.tsv.gz", sep="\t", index=False, compression="gzip")

cols = ["chr", "start", "end", "var_id", "g_chr", "g_start", "g_end", "gene_id", "dist"]
reg = pd.concat([pd.read_csv(a.reg_snps, sep="\t", header=None, names=cols),
                 pd.read_csv(a.reg_indels, sep="\t", header=None, names=cols)],
                ignore_index=True).drop_duplicates(subset=["chr", "start", "var_id"])
cand = pd.read_csv(a.candidates, sep="\t")
cv = reg[reg.gene_id.isin(cand.gene_id)].copy()
# var_id is CHROM:POS:REF:ALT with a chr prefix; the genotype table has no prefix
cv["key"] = [tuple([p[0].replace("chr", "", 1)] + p[1:]) for p in cv.var_id.str.split(":")]
idx = allv.set_index(["chrom", "pos", "ref", "alt"])
hit = cv[cv.key.isin(idx.index)]
k_c = int(idx.loc[list(hit.key), "concordant"].sum()) if len(hit) else 0

n, k = len(allv), int(allv.concordant.sum())
p0 = k / n
P = binomtest(k_c, len(hit), p0, alternative="greater").pvalue if len(hit) else float("nan")
summary = pd.DataFrame([dict(
    definition="kidney-mass ranking, C57BL/6J fixed at 0",
    background_n=n, background_concordant=k, background_pct=round(100 * p0, 2),
    candidate_genes=len(cand), candidate_variants_in_set=len(hit),
    candidate_variants_total=len(cv), candidate_concordant=k_c,
    candidate_pct=round(100 * k_c / len(hit), 1) if len(hit) else None,
    binom_one_sided_P=round(P, 4))])
summary.to_csv(f"{a.outdir}/summary.tsv", sep="\t", index=False)
pd.DataFrame(audit).to_csv(f"{a.outdir}/input_audit.tsv", sep="\t", index=False)
print(pd.DataFrame(audit).to_string(index=False)); print()
print(summary.to_string(index=False))
