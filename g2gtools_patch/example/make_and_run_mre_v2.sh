#!/bin/bash
# Item 4 minimal reproducible example for the g2gtools 0.2.7 BAM reverse-conversion defect.
# Builds a 3-kb synthetic reference with one SNP, one 3-bp deletion and one 3-bp insertion
# (strain "FVB_NJ", GT 1/1), constructs the personalised genome with the same g2gtools
# sub-commands as run_phase2_condB.sh, places reads on personalised-genome coordinates,
# and converts them back to reference coordinates with (a) the installed g2gtools CLI and
# (b) the patched wrapper g2g_bam_convert.py + bsam_fixed.py used in run_quantify_cleanup.sh.
set -uo pipefail
cd "$(dirname "$0")"
rm -rf work && mkdir work && cd work
BCF=/usr/bin/bcftools; SAMT=/usr/local/bin/samtools; PY=/usr/bin/python3
$PY - <<'P'
import random
random.seed(1)
seq = "".join(random.choice("ACGT") for _ in range(3000))
seq2 = "".join(random.choice("ACGT") for _ in range(1000))  # chr2: contig WITHOUT any variant (like chrM/chrY/scaffolds for FVB_NJ)
open("ref.fa","w").write(">chr1\n" + "\n".join(seq[i:i+60] for i in range(0,len(seq),60)) + "\n>chr2\n" + "\n".join(seq2[i:i+60] for i in range(0,len(seq2),60)) + "\n")
# variants (1-based VCF POS)
snp_pos = 800;  snp_ref = seq[snp_pos-1]; snp_alt = {"A":"G","C":"T","G":"A","T":"C"}[snp_ref]
del_pos = 1500; del_ref = seq[del_pos-1:del_pos-1+4]; del_alt = seq[del_pos-1]            # deletes 3 bp after POS
ins_pos = 2200; ins_ref = seq[ins_pos-1]; ins_alt = ins_ref + "TTT"                        # inserts 3 bp after POS
with open("fvb.vcf","w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=3000>\n##contig=<ID=chr2,length=1000>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFVB_NJ\n")
    for p,r,a in [(snp_pos,snp_ref,snp_alt),(del_pos,del_ref,del_alt),(ins_pos,ins_ref,ins_alt)]:
        f.write(f"chr1\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t1/1\n")
P
$BCF view fvb.vcf -Oz -o fvb.vcf.gz && $BCF index -t fvb.vcf.gz && $SAMT faidx ref.fa
g2gtools vcf2vci -f ref.fa -i fvb.vcf.gz -s FVB_NJ -o fvb.vci > vcf2vci.log 2>&1
g2gtools patch -i ref.fa -c fvb.vci.gz -o fvb.patched.fa > patch.log 2>&1
g2gtools transform -i fvb.patched.fa -c fvb.vci.gz -o fvb.fa > transform.log 2>&1
$SAMT faidx fvb.fa
cat fvb.fa.fai
# reads placed on FVB coordinates (0-based starts): before deletion, between deletion and insertion, after insertion
$PY - <<'P'
import pysam
fvb = pysam.FastaFile("fvb.fa"); L = fvb.get_reference_length("chr1")
hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": L}, {"SN": "chr2", "LN": fvb.get_reference_length("chr2")}]}
starts = [100, 750, 1700, 2600]  # FVB-coordinate 0-based starts (750 overlaps the SNP at ref POS 800)
with pysam.AlignmentFile("fvb_coords.sam", "w", header=hdr) as out:
    for i, s in enumerate(starts):
        a = pysam.AlignedSegment(); a.query_name = f"read{i}"; a.flag = 0
        a.reference_id = 0; a.reference_start = s; a.mapping_quality = 60
        a.cigar = [(0, 100)]; a.query_sequence = fvb.fetch("chr1", s, s + 100); a.query_qualities = pysam.qualitystring_to_array("I" * 100)
        out.write(a)
    b = pysam.AlignedSegment(); b.query_name = "read_chr2"; b.flag = 0; b.reference_id = 1; b.reference_start = 200; b.mapping_quality = 60
    b.cigar = [(0, 100)]; b.query_sequence = fvb.fetch("chr2", 200, 300); b.query_qualities = pysam.qualitystring_to_array("I" * 100); out.write(b)
open("expected.tsv","w").write("read\tfvb_start0\texpected_ref_start0\n" + "\n".join(
    f"read{i}\t{s}\t{s + (3 if 1497 <= s < 2200 else 0)}" for i, s in enumerate(starts)) + "\nread_chr2\t200\t200\n")
P
$SAMT sort -o fvb_coords.bam fvb_coords.sam && $SAMT index fvb_coords.bam
echo "================ (a) installed g2gtools CLI: g2gtools convert -f bam --reverse ================"
g2gtools convert -i fvb_coords.bam -c fvb.vci.gz -f bam --reverse -o cli_reverted.bam > cli_convert.log 2>&1; echo "exit code: $?"
tail -25 cli_convert.log
ls -la cli_reverted.bam* 2>/dev/null
[ -s cli_reverted.bam ] && { $SAMT view -H cli_reverted.bam | grep '^@SQ'; $SAMT view cli_reverted.bam | cut -f1-6; }
echo "================ (a2) installed bsam.convert_bam_file with a VCIFile that is constructed but not parsed ================"
$PY - > unparsed_vci.log 2>&1 <<'P'
from g2gtools import vci, bsam
v = vci.VCIFile("fvb.vci.gz", reverse=True)
print("mapping_tree size before parse:", len(v.mapping_tree))
try:
    bsam.convert_bam_file(v, "fvb_coords.bam", "installed_unparsed.bam", reverse=True)
except Exception as e:
    print("EXCEPTION:", type(e).__name__, e)
P
echo "exit: $?"; tail -8 unparsed_vci.log
echo "================ (a3) installed bsam.convert_bam_file with a PARSED VCIFile (isolates the SQ-header defect) ================"
$PY - > installed_parsed_vci.log 2>&1 <<'P'
from g2gtools import vci, bsam
v = vci.VCIFile("fvb.vci.gz", reverse=True); v.parse(True)
print("mapping_tree size after parse:", len(v.mapping_tree))
try:
    bsam.convert_bam_file(v, "fvb_coords.bam", "installed_parsed.bam", reverse=True)
except Exception as e:
    import traceback; traceback.print_exc(); print("EXCEPTION:", type(e).__name__, e)
P
tail -12 installed_parsed_vci.log
echo "================ (b) patched: g2g_bam_convert.py + bsam_fixed.py ================"
cp ../../g2g_bam_convert_v2.py ../../bsam_fixed.py .
$PY g2g_bam_convert_v2.py -c fvb.vci.gz -i fvb_coords.bam -o patched_reverted.bam --reverse > patched_convert.log 2>&1; echo "exit code: $?"
cat patched_convert.log | tail -5
$SAMT view -H patched_reverted.bam | grep '^@SQ'
$SAMT view patched_reverted.bam | cut -f1-6 > patched_reverted.tsv; cat patched_reverted.tsv
echo "---- expected ----"; cat expected.tsv
$PY - <<'P'
import pysam
exp = {l.split()[0]: int(l.split()[2]) for l in open("expected.tsv").read().splitlines()[1:]}
print("reads written to patched_reverted.bam.unmapped:", [x.query_name for x in pysam.AlignmentFile("patched_reverted.bam.unmapped", check_sq=False)])
ref = pysam.FastaFile("ref.fa"); ok = True
for a in pysam.AlignmentFile("patched_reverted.bam"):
    mism = sum(x != y for x, y in zip(a.query_sequence, ref.fetch(a.reference_name, a.reference_start, a.reference_start + 100)))
    good = a.reference_start == exp[a.query_name] and mism <= 1
    print(a.query_name, a.reference_name, "ref_start0", a.reference_start, "expected", exp[a.query_name], "cigar", a.cigarstring, "mismatches vs reference (1 expected only for the SNP read):", mism, "->", "OK" if good else "CHECK")
    ok &= good
print("ALL_OK" if ok else "NOT_ALL_OK")
P
