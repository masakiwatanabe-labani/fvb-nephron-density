#!/bin/bash
set -euo pipefail
BCFTOOLS=/usr/bin/bcftools
cd /usr/local/jupyter/FVB_B6_glo
for v in snps indels; do
  echo "[$(date +%T)] norm: $v"
  $BCFTOOLS norm -f refs/GRCm39.primary_assembly.genome.fa -m -any \
    results/phase1/variants_chr/FVB.${v}.chr.vcf.gz -Oz -o results/phase2/g2g/FVB.${v}.norm.vcf.gz
  $BCFTOOLS index -t results/phase2/g2g/FVB.${v}.norm.vcf.gz
  echo "[$(date +%T)] done: $v ($($BCFTOOLS view -H results/phase2/g2g/FVB.${v}.norm.vcf.gz | wc -l) records)"
done
echo "ALL_NORM_DONE"
