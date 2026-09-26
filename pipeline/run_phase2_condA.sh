#!/bin/bash
set -euo pipefail
cd /usr/local/jupyter/FVB_B6_glo
echo "[$(date +%T)] vcf2vci (condition A: SNP-only)"
g2gtools vcf2vci -f refs/GRCm39.primary_assembly.genome.fa \
  -i results/phase2/g2g/FVB.snps.norm.vcf.gz \
  -s FVB_NJ -p 8 \
  -o results/phase2/g2g/condA/FVB_snponly.vci
echo "[$(date +%T)] patch"
g2gtools patch -i refs/GRCm39.primary_assembly.genome.fa \
  -c results/phase2/g2g/condA/FVB_snponly.vci.gz \
  -p 8 \
  -o results/phase2/g2g/condA/FVB_snponly.fa
cp refs/gencode.vM39.annotation.gtf results/phase2/g2g/condA/FVB_snponly.gtf
echo "ALL_CONDA_DONE"
