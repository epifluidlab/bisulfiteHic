#!/usr/bin/env python3
"""Verify SEQ / QUAL / strand consistency of a Bhmem BAM against its input FASTQ
and the reference. Reads SAM text on stdin (samtools view), standard library only.

Usage:
  samtools view out.bam | check_bam.py --ref genome.fa --r1 R1.fq --r2 R2.fq \
      [--pbat] [--min-pairs N]

Checks (each prints PASS/FAIL; exit status 1 if any fails):
  qual_orientation   QUAL is reversed exactly when SEQ is the reverse complement
                     of the FASTQ read.
  seq_vs_reference   Aligned bases agree with the reference apart from bisulfite
                     conversions (ref C/read T, ref G/read A).
  strand_flag        The 0x10 flag agrees with the SEQ orientation (inverted for
                     PBAT, where Bhmem flips the flag but keeps SEQ on the
                     reference forward strand).
  pairs_written      At least --min-pairs read pairs are present.
"""
import argparse
import gzip
import re
import sys

COMP = str.maketrans("ACGTN", "TGCAN")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def read_fastq(path):
    out = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            head = fh.readline()
            if not head:
                break
            seq = fh.readline().strip()
            fh.readline()
            qual = fh.readline().strip()
            out[head[1:].split()[0]] = (seq, qual)
    return out


def non_bisulfite_mismatch_fraction(seq, cigar, chrom_seq, pos):
    """Fraction of aligned (M) bases that differ from the reference in a way a
    bisulfite conversion cannot explain. `pos` is 1-based."""
    qi, ri, aligned, bad = 0, pos - 1, 0, 0
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in "M=X":
            for k in range(n):
                r, q = chrom_seq[ri + k], seq[qi + k]
                aligned += 1
                if q != r and not ((r == "C" and q == "T") or (r == "G" and q == "A")):
                    bad += 1
            qi += n
            ri += n
        elif op in "IS":
            qi += n
        elif op in "DN":
            ri += n
    return (bad / aligned) if aligned else 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--pbat", action="store_true")
    ap.add_argument("--min-pairs", type=int, default=1)
    ap.add_argument("--max-mismatch-frac", type=float, default=0.15)
    args = ap.parse_args()

    ref = read_fasta(args.ref)
    fq = {1: read_fastq(args.r1), 2: read_fastq(args.r2)}

    n = 0
    names = set()
    fails = {"qual_orientation": [], "seq_vs_reference": [], "strand_flag": []}
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        name, flag, chrom, pos, cigar, seq, qual = f[0], int(f[1]), f[2], int(f[3]), f[5], f[9], f[10]
        end = 1 if flag & 0x40 else 2
        fseq, fqual = fq[end][name]
        n += 1
        names.add(name)
        tag = f"{name}/{end} flag={flag} {chrom}:{pos} {cigar}"
        if "H" in cigar:
            # Hard-clipped records cannot be compared base-by-base with the FASTQ.
            continue

        is_fwd, is_rc = seq == fseq, seq == rc(fseq)
        if not (is_fwd or is_rc):
            fails["qual_orientation"].append(tag + " SEQ is neither the read nor its reverse complement")
            continue
        if is_fwd and is_rc:
            continue  # reverse-palindromic read: orientation undefined
        expect_qual = fqual if is_fwd else fqual[::-1]
        if qual != expect_qual:
            fails["qual_orientation"].append(
                tag + f" SEQ is {'forward' if is_fwd else 'revcomp'} but QUAL is not in the matching order")

        frac = non_bisulfite_mismatch_fraction(seq, cigar, ref[chrom], pos)
        if frac > args.max_mismatch_frac:
            fails["seq_vs_reference"].append(tag + f" non-bisulfite mismatch fraction {frac:.2f}")

        flag_rev = bool(flag & 0x10)
        expect_rev = (not is_rc) if args.pbat else is_rc
        if flag_rev != expect_rev:
            fails["strand_flag"].append(tag + f" 0x10={int(flag_rev)} but SEQ is {'revcomp' if is_rc else 'forward'}")

    ok = True
    for check, bad in fails.items():
        status = "PASS" if not bad else "FAIL"
        ok &= not bad
        print(f"  {status}  {check}: {len(bad)} of {n} records inconsistent")
        for msg in bad[:3]:
            print(f"        e.g. {msg}")
    enough = len(names) >= args.min_pairs
    ok &= enough
    print(f"  {'PASS' if enough else 'FAIL'}  pairs_written: {len(names)} pairs (need >= {args.min_pairs})")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
