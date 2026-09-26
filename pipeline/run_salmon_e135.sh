#!/bin/bash
set -uo pipefail
source /root/miniforge3/etc/profile.d/conda.sh
conda activate fvb_pipeline
cd /usr/local/jupyter/FVB_B6_glo
mkdir -p results/phase0/sex_check
for sample in B6135-1 B6135-2 B6135-3 FVB135-1 FVB135-2 FVB135-3; do
  echo "[$(date +%T)] salmon quant: $sample"
  salmon quant -i refs/salmon_index -l A \
    -1 "results/phase0/fastp/${sample}_trimmed_R1.fastq.gz" \
    -2 "results/phase0/fastp/${sample}_trimmed_R2.fastq.gz" \
    -p 4 --validateMappings \
    -o "results/phase0/sex_check/${sample}" \
    > "results/phase0/sex_check/${sample}.log" 2>&1
done
echo "ALL_SALMON_E135_DONE"
