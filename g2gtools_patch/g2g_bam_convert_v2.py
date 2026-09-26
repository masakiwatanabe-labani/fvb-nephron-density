#!/usr/bin/env python3
"""Reverse (personalised -> reference) BAM coordinate conversion with g2gtools 0.2.7.

v2 of g2g_bam_convert.py. Differences from v1:
  * contigs that are declared in the VCI header (##CONTIG) but have no VCI records are added to the
    mapping tree as a single identity interval (0, length) -> (0, length). Without this, g2gtools
    leaves them out of the tree and every alignment on them (for FVB_NJ: chrM, chrY and the unplaced
    scaffolds) is written to the .unmapped file.
  * the list of contigs given an identity mapping is written to stderr.
bsam_fixed.py (rebuilds @SQ from the VCI ##CONTIG lengths, closes files) is unchanged. As in v1 the
conversion is only valid for --reverse, because @SQ lengths are the reference lengths."""
import argparse, os, sys
from bx.intervals.intersection import Interval, IntervalTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g2gtools import vci
import bsam_fixed as bsam

def add_identity_contigs(vci_file):
    added = []
    for contig, length in vci_file.contigs.items():
        if contig not in vci_file.mapping_tree:
            tree = IntervalTree()
            tree.insert_interval(Interval(0, length, vci.IntervalInfo(contig, 0, length, None, None, None, None)))
            vci_file.mapping_tree[contig] = tree
            added.append(contig)
    return added

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--vci", required=True); ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True); ap.add_argument("--reverse", action="store_true")
    a = ap.parse_args()
    if not a.reverse:
        sys.exit("only --reverse is supported (the @SQ header is built from reference contig lengths)")
    v = vci.VCIFile(a.vci, reverse=True)
    v.parse(True)
    n_var = len(v.mapping_tree)
    added = add_identity_contigs(v)
    print(f"Parsed VCI: {n_var} contigs with variant records; identity mapping added for {len(added)} contigs: {','.join(added)}", file=sys.stderr)
    bsam.convert_bam_file(v, a.input, a.output, reverse=True)

if __name__ == "__main__":
    main()
