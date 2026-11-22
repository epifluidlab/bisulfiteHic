#!/usr/bin/env python3
"""Verify that Bhmem's -enzymeList changes which candidate pair is chosen.

The `enz` test reads have an R2 that splits 50/50 between locus A and locus B
(both on chrA), with tied mapping quality. Run once with an enzyme list that
supports locus A and once with a list that supports locus B; the split read must
be placed at the supported locus in each run.

Usage: check_enzyme.py enz_truth.tsv runA.sam runB.sam [--min-frac 0.9]
"""
import sys


def split_read_positions(sam_path):
    """name -> 1-based POS of the second-in-pair record (the split read)."""
    pos = {}
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t")
            if int(f[1]) & 0x80:
                pos[f[0]] = (f[2], int(f[3]))
    return pos


def near(pos, chrom, start, span=60):
    return pos is not None and pos[0] == chrom and start - span <= pos[1] - 1 <= start + span


def main():
    truth_path, sam_a, sam_b = sys.argv[1:4]
    min_frac = float(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[4] == "--min-frac" else 0.9
    pos_a, pos_b = split_read_positions(sam_a), split_read_positions(sam_b)
    total = follow = 0
    examples = []
    with open(truth_path) as fh:
        for line in fh:
            name, chrom, sa, sb = line.split()
            sa, sb = int(sa), int(sb)
            total += 1
            at_a = near(pos_a.get(name), chrom, sa)
            at_b = near(pos_b.get(name), chrom, sb)
            if at_a and at_b:
                follow += 1
            elif len(examples) < 3:
                examples.append(f"{name}: locusA={sa} locusB={sb} listA->{pos_a.get(name)} listB->{pos_b.get(name)}")
    frac = follow / total if total else 0.0
    ok = frac >= min_frac
    print(f"  {'PASS' if ok else 'FAIL'}  enzyme_list_respected: {follow} of {total} split reads "
          f"follow the enzyme list in both runs (need >= {min_frac:.0%})")
    for e in examples:
        print(f"        e.g. {e}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
