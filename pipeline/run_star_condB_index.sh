#!/bin/bash
set -euo pipefail
cd /usr/local/jupyter/FVB_B6_glo
STAR --runMode genomeGenerate \
  --runThreadN 8 \
  --genomeDir refs/star_condB \
  --genomeFastaFiles results/phase2/g2g/condB/FVB_full.fa \
  --sjdbGTFfile results/phase2/g2g/condB/FVB_full.gtf \
  --sjdbOverhang 149 \
  --genomeSAindexNbases 14
echo "STAR_CONDB_INDEX_DONE"
