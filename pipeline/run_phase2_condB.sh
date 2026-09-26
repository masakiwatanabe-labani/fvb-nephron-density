#!/bin/bash
set -euo pipefail
cd /usr/local/jupyter/FVB_B6_glo
echo "[$(date +%T)] vcf2vci (condition B: SNP+indel)"
g2gtools vcf2vci -f refs/GRCm39.primary_assembly.genome.fa \
  -i results/phase2/g2g/FVB.snps.norm.vcf.gz \
  -i results/phase2/g2g/FVB.indels.norm.vcf.gz \
  -s FVB_NJ -p 8 \
  -o results/phase2/g2g/condB/FVB_full.vci
echo "[$(date +%T)] patch"
g2gtools patch -i refs/GRCm39.primary_assembly.genome.fa \
  -c results/phase2/g2g/condB/FVB_full.vci.gz \
  -p 8 \
  -o results/phase2/g2g/condB/FVB_full.patched.fa
echo "[$(date +%T)] transform"
g2gtools transform -i results/phase2/g2g/condB/FVB_full.patched.fa \
  -c results/phase2/g2g/condB/FVB_full.vci.gz \
  -p 8 \
  -o results/phase2/g2g/condB/FVB_full.fa
echo "[$(date +%T)] convert (GTF)"
g2gtools convert -i refs/gencode.vM39.annotation.gtf \
  -c results/phase2/g2g/condB/FVB_full.vci.gz \
  -f gtf \
  -o results/phase2/g2g/condB/FVB_full.gtf
echo "ALL_CONDB_DONE"
