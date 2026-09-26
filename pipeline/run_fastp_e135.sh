#!/bin/bash
set -uo pipefail
source /root/miniforge3/etc/profile.d/conda.sh
conda activate fvb_pipeline
cd /usr/local/jupyter/FVB_B6_glo
for sample in B6135-1 B6135-2 B6135-3 FVB135-1 FVB135-2 FVB135-3; do
  echo "[$(date +%T)] fastp: $sample"
  r1="input/${sample}_S1_L001_R1_001.fastq.gz"
  r2="input/${sample}_S1_L001_R2_001.fastq.gz"
  fastp \
    -i "$r1" -I "$r2" \
    -o "results/phase0/fastp/${sample}_trimmed_R1.fastq.gz" \
    -O "results/phase0/fastp/${sample}_trimmed_R2.fastq.gz" \
    -j "results/phase0/fastp/${sample}.fastp.json" \
    -h "results/phase0/fastp/${sample}.fastp.html" \
    -w 1 \
    2> "results/phase0/fastp/${sample}.fastp.log"
done
echo "ALL_FASTP_E135_DONE"
