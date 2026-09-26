#!/bin/bash
# Four-strain genotypes for exactly the records of the allele-wise private-variant set.
#
# Defect in rerun3/F7/extract_concordance.sh that this fixes:
#   `bcftools view -R <private.vcf.gz>` treats the private VCF as a region list, and for an indel a
#   region spans more than one base, so every source record overlapping that span is returned -
#   including records that are not private variants. For the indel file this turned 10,912 private
#   records into 28,462 rows covering 14,045 positions, 3,133 of which are not private variants.
# Here the source records are decomposed to one alternate allele each, exactly as the private set was
# built, so that CHROM/POS/REF/ALT identifies a record and the join is exact.
#
# Usage: extract_genotypes_v3.sh <priv_snps.vcf.gz> <priv_indels.vcf.gz> <mgp_dir> <outdir>
set -euo pipefail
PSNP=$1; PIND=$2; MGP=$3; OUT=$4
mkdir -p "$OUT"
for v in snps indels; do
  src=$([ "$v" = snps ] && echo "$PSNP" || echo "$PIND")
  echo "[$(date +%T)] $v"
  bcftools view -R "$src" -s BALB_cJ,DBA_2J,C57BL_6NJ,FVB_NJ "$MGP/mgp_REL2021_${v}.vcf.gz" -Ou 2>/dev/null \
    | bcftools norm -m -any -Ou 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT,]\n' > "$OUT/${v}_genotypes_v3.txt"
  echo "[$(date +%T)] $v done: $(wc -l < "$OUT/${v}_genotypes_v3.txt") rows"
done
echo EXTRACT_V3_DONE
