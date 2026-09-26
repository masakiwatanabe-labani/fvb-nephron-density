#!/bin/bash
# Re-count an arm with the library-appropriate strandedness and fragment-level (pair) counting.
#
#   -s 2              reverse-stranded. salmon `-l A` auto-detected ISR for all 12 libraries
#                     (results/phase0/sex_check/<sample>/lib_format_counts.json, compatible_fragment_ratio 1.0).
#   --countReadPairs  count fragments, not the two reads of a pair separately. `-p` alone declares the
#                     input paired-end; in subread 2.x it does not by itself switch counting to fragments.
#
# Usage: featurecounts_s2_pairs.sh <project_root> <bam_pattern> <gtf> <outdir> <tag> [parallel_jobs]
#   <bam_pattern> contains SAMPLE where the sample name goes, e.g.
#     results/phase2/star_uncorrected/SAMPLE/Aligned.sortedByCoord.out.bam
set -euo pipefail
ROOT=$1; BAMPAT=$2; GTF=$3; OUT=$4; TAG=$5; J=${6:-2}
FC=${FEATURECOUNTS:-featureCounts}
mkdir -p "$OUT"
$FC -v 2>&1 | grep -o 'v[0-9.]*' | head -1 > "$OUT/featureCounts_version.txt" || true

run_one() {
  s=$1
  bam=${BAMPAT//SAMPLE/$s}
  cd "$ROOT"
  "$FC" -T 4 -p --countReadPairs -s 2 -a "$GTF" -t exon -g gene_id \
    -o "$OUT/$s.$TAG.counts.txt" "$bam" > "$OUT/$s.$TAG.featureCounts.log" 2>&1
  echo "done $s"
}
export -f run_one
export ROOT BAMPAT GTF OUT TAG FC

samples="B6135-1 B6135-2 B6135-3 B6P1-1 B6P1-2 B6P1-3 FVB135-1 FVB135-2 FVB135-3 FVBP1-1 FVBP1-2 FVBP1-3"
printf "%s\n" $samples | xargs -P "$J" -I{} bash -c 'run_one "$@"' _ {}
echo "ALL_DONE $TAG"
