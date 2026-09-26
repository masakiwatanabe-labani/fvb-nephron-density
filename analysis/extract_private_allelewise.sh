#!/bin/bash
# Allele-level extraction of FVB/NJ private variants.
#
# Defect in run_phase1_variants_v2.sh that this fixes:
#   the record was kept if  GT[FVB]="alt"  AND  INFO/AC<=2 .
#   On a multi-allelic record INFO/AC is a vector (one count per ALT) and bcftools evaluates the
#   comparison as true when ANY element satisfies it. The two conditions were therefore not tied to
#   the same allele: a site passed when some OTHER alt allele was rare, whatever allele FVB carried.
#   Decomposing to biallelic records before computing AC ties them together.
#
# Usage: extract_private_allelewise.sh <snps.vcf.gz> <indels.vcf.gz> <outdir> [fvb_index] [max_ac]
set -euo pipefail
SNPS="$1"; INDELS="$2"; OUT="$3"; IDX="${4:-24}"; MAXAC="${5:-2}"
mkdir -p "$OUT"

# FILTER=PASS is meaningful for the SNP file; the indel file has FILTER uniformly "." in this release.
echo "[$(date +%T)] SNPs: decompose -> per-allele AC over all strains -> FVB carries that allele AND AC<=${MAXAC}"
bcftools view -f PASS "$SNPS" -Ou \
  | bcftools norm -m -any -Ou \
  | bcftools +fill-tags -- -t AC \
  | bcftools view -i "GT[${IDX}]=\"alt\" && INFO/AC<=${MAXAC}" -Ou \
  | bcftools view -s FVB_NJ -Oz -o "$OUT/private_snps_allelewise.vcf.gz"
bcftools index -t "$OUT/private_snps_allelewise.vcf.gz"
echo "[$(date +%T)] SNPs done: $(bcftools view -H "$OUT/private_snps_allelewise.vcf.gz" | wc -l)"

echo "[$(date +%T)] indels"
bcftools view "$INDELS" -Ou \
  | bcftools norm -m -any -Ou \
  | bcftools +fill-tags -- -t AC \
  | bcftools view -i "GT[${IDX}]=\"alt\" && INFO/AC<=${MAXAC}" -Ou \
  | bcftools view -s FVB_NJ -Oz -o "$OUT/private_indels_allelewise.vcf.gz"
bcftools index -t "$OUT/private_indels_allelewise.vcf.gz"
echo "[$(date +%T)] indels done: $(bcftools view -H "$OUT/private_indels_allelewise.vcf.gz" | wc -l)"
echo ALLELEWISE_EXTRACTION_DONE
