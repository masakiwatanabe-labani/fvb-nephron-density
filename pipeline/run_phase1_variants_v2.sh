#!/bin/bash
set -euo pipefail
BCFTOOLS=/usr/bin/bcftools
REFS=/usr/local/jupyter/FVB_B6_glo/refs
OUT=/usr/local/jupyter/FVB_B6_glo/results/phase1/variants
cd "$OUT"

FVB_IDX_SNPS=24
FVB_IDX_INDELS=24

echo "[$(date +%T)] Step A: FVB_NJ genome-wide extraction (for Phase 2 personal genome)"
# SNPs: FILTER=PASS is meaningful here
$BCFTOOLS view -s FVB_NJ -f PASS -c1 "$REFS/mgp_REL2021_snps.vcf.gz" -Oz -o "FVB.snps.vcf.gz"
$BCFTOOLS index -t "FVB.snps.vcf.gz"
echo "[$(date +%T)] FVB.snps.vcf.gz done: $($BCFTOOLS view -H FVB.snps.vcf.gz | wc -l) variants"

# indels: FILTER is uniformly "." in this MGP release (verified: no PASS/LowQual anywhere) -> drop -f
$BCFTOOLS view -s FVB_NJ -c1 "$REFS/mgp_REL2021_indels.vcf.gz" -Oz -o "FVB.indels.vcf.gz"
$BCFTOOLS index -t "FVB.indels.vcf.gz"
echo "[$(date +%T)] FVB.indels.vcf.gz done: $($BCFTOOLS view -H FVB.indels.vcf.gz | wc -l) variants"

echo "[$(date +%T)] Step B: genome-wide private variant extraction"
echo "  (population AC computed on all 52 strains BEFORE any sample subsetting;"
echo "   site kept only if FVB_NJ itself carries the ALT allele AND population AC<=2;"
echo "   -s applied only as the final column-trim, after filtering, so it cannot"
echo "   silently recompute/overwrite AC)"

# SNPs
$BCFTOOLS view -f PASS "$REFS/mgp_REL2021_snps.vcf.gz" -Ou \
  | $BCFTOOLS +fill-tags -- -t AC \
  | $BCFTOOLS view -i "GT[${FVB_IDX_SNPS}]=\"alt\" && INFO/AC<=2" -Ou \
  | $BCFTOOLS view -s FVB_NJ -Oz -o "private_snps.vcf.gz"
$BCFTOOLS index -t "private_snps.vcf.gz"
echo "[$(date +%T)] private_snps.vcf.gz done: $($BCFTOOLS view -H private_snps.vcf.gz | wc -l) variants"

# indels (no -f PASS: not meaningful for this file)
$BCFTOOLS view "$REFS/mgp_REL2021_indels.vcf.gz" -Ou \
  | $BCFTOOLS +fill-tags -- -t AC \
  | $BCFTOOLS view -i "GT[${FVB_IDX_INDELS}]=\"alt\" && INFO/AC<=2" -Ou \
  | $BCFTOOLS view -s FVB_NJ -Oz -o "private_indels.vcf.gz"
$BCFTOOLS index -t "private_indels.vcf.gz"
echo "[$(date +%T)] private_indels.vcf.gz done: $($BCFTOOLS view -H private_indels.vcf.gz | wc -l) variants"

echo "ALL_PHASE1_VARIANTS_V2_DONE"
