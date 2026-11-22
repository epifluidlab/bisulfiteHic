#!/usr/bin/env python3
"""Generate a tiny synthetic reference and bisulfite Hi-C read pairs for the
Bhmem regression tests. Pure standard library, fully seeded (reproducible).

Outputs under OUTDIR:
  ref/genome.fa                       two random chromosomes (chrA, chrB)
  ref/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa
  ref/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa
  reads/dir_R{1,2}.fq                 directional library (R1 = C->T strand)
  reads/pbat_R{1,2}.fq                PBAT library (R1 = complement strand)
  reads/enz_R{1,2}.fq                 PBAT pairs whose R1 is a 50/50 chimera
  reads/enz_truth.tsv                 read name, chrom, locusA start, locusB start
  reads/enz_A.bed, reads/enz_B.bed    enzyme-site lists favouring locus A / B

Every read gets a strictly non-palindromic quality string (high at the 5' end,
low at the 3' end) so the orientation of QUAL in the BAM can be verified.
"""
import os
import random
import sys

READ_LEN = 100
CHROMS = {"chrA": 60000, "chrB": 60000}
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def bisulfite(seq, rng, meth=0.3):
    """C->T conversion of one strand; each C is 'methylated' (kept) with prob meth."""
    return "".join(("C" if rng.random() < meth else "T") if b == "C" else b for b in seq)


def qual_string(n):
    """Monotonically decreasing Phred string: 'J' (Q41) at 5' end down to '#' (Q2)."""
    return "".join(chr(33 + 41 - (i * 39) // max(n - 1, 1)) for i in range(n))


def write_fasta(path, records, width=60):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def reads_from_molecule(mol, top, pbat, rng):
    """Return (R1, R2) for a double-stranded molecule whose top strand is `mol`.

    One strand is bisulfite converted (top or bottom). Directional libraries read
    the converted strand as R1; PBAT libraries read its complement as R1.
    """
    bs = bisulfite(mol if top else rc(mol), rng)
    n = min(READ_LEN, len(bs))
    if pbat:
        return rc(bs)[:n], bs[:n]
    return bs[:n], rc(bs)[:n]


def segment(genome, chrom, start, length, strand):
    seq = genome[chrom][start:start + length]
    return seq if strand == "+" else rc(seq)


def random_locus(rng, length, chrom=None, avoid=None):
    """Pick (chrom, start); if `avoid`=(chrom,start) keep >=5 kb away on the same chrom."""
    while True:
        c = chrom or rng.choice(sorted(CHROMS))
        s = rng.randrange(500, CHROMS[c] - length - 500)
        if avoid and avoid[0] == c and abs(avoid[1] - s) < 5000:
            continue
        return c, s


def make_library(genome, rng, pbat, n_mol):
    """Mixed library: plain paired-end fragments plus Hi-C ligation chimeras."""
    kinds = [
        ("pe", 250, 0),        # single locus, no ligation
        ("clean", 150, 150),   # junction outside both reads
        ("junc5", 60, 190),    # junction inside the read from the first segment
        ("junc3", 190, 60),    # junction inside the read from the second segment
        ("both", 70, 70),      # short molecule, junction inside both reads
    ]
    pairs = []
    for i in range(n_mol):
        kind, len_a, len_b = kinds[i % len(kinds)]
        ca, sa = random_locus(rng, len_a)
        mol = segment(genome, ca, sa, len_a, rng.choice("+-"))
        if len_b:
            same_chrom = rng.random() < 0.5
            cb, sb = random_locus(rng, len_b, chrom=ca if same_chrom else None, avoid=(ca, sa))
            mol += segment(genome, cb, sb, len_b, rng.choice("+-"))
        r1, r2 = reads_from_molecule(mol, rng.random() < 0.5, pbat, rng)
        pairs.append((f"{'pbat' if pbat else 'dir'}_{kind}_{i:05d}", r1, r2))
    return pairs


def make_enzyme_set(genome, rng, n_mol):
    """PBAT pairs designed so two candidate pairs tie on MAPQ.

    Molecule = A[0:50] + B[0:250], both loci on chrA, forward strand, top strand
    converted. R2 (bs[:100]) = 50 bp of A then 50 bp of B  -> split alignment.
    R1 (rc(bs)[:100])        = B[150:250]                   -> single alignment.
    Enzyme list "A" has sites at A+20 and B+200 (supports placing the split read
    at A); list "B" has sites at B+20 and B+200 (supports placing it at B).
    """
    pairs, truth, bed_a, bed_b = [], [], [], []
    # Loci come from a shuffled 1 kb grid so no two molecules' sites can interfere
    # (and the generator cannot loop forever looking for free space).
    grid = list(range(2000, CHROMS["chrA"] - 2000, 1000))
    rng.shuffle(grid)
    if 2 * n_mol > len(grid):
        raise ValueError("too many enzyme-test molecules for the chromosome size")
    for i in range(n_mol):
        sa, sb = grid[2 * i], grid[2 * i + 1]
        mol = genome["chrA"][sa:sa + 50] + genome["chrA"][sb:sb + 250]
        r1, r2 = reads_from_molecule(mol, True, True, rng)
        name = f"enz_{i:05d}"
        pairs.append((name, r1, r2))
        truth.append((name, "chrA", sa, sb))
        bed_a += [("chrA", sa + 20, sa + 24), ("chrA", sb + 200, sb + 204)]
        bed_b += [("chrA", sb + 20, sb + 24), ("chrA", sb + 200, sb + 204)]
    return pairs, truth, sorted(bed_a), sorted(bed_b)


def write_fastq(prefix, pairs):
    with open(prefix + "_R1.fq", "w") as f1, open(prefix + "_R2.fq", "w") as f2:
        for name, r1, r2 in pairs:
            f1.write(f"@{name}\n{r1}\n+\n{qual_string(len(r1))}\n")
            f2.write(f"@{name}\n{r2}\n+\n{qual_string(len(r2))}\n")


def main():
    outdir = sys.argv[1]
    rng = random.Random(20261121)
    genome = {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in CHROMS.items()}

    ref = os.path.join(outdir, "ref")
    ct = os.path.join(ref, "Bisulfite_Genome", "CT_conversion")
    ga = os.path.join(ref, "Bisulfite_Genome", "GA_conversion")
    reads = os.path.join(outdir, "reads")
    for d in (ct, ga, reads):
        os.makedirs(d, exist_ok=True)

    order = sorted(genome)
    write_fasta(os.path.join(ref, "genome.fa"), [(c, genome[c]) for c in order])
    # Same naming convention as bismark_genome_preparation.
    write_fasta(os.path.join(ct, "genome_mfa.CT_conversion.fa"),
                [(c + "_CT_converted", genome[c].replace("C", "T")) for c in order])
    write_fasta(os.path.join(ga, "genome_mfa.GA_conversion.fa"),
                [(c + "_GA_converted", genome[c].replace("G", "A")) for c in order])

    write_fastq(os.path.join(reads, "dir"), make_library(genome, rng, False, 500))
    write_fastq(os.path.join(reads, "pbat"), make_library(genome, rng, True, 500))

    pairs, truth, bed_a, bed_b = make_enzyme_set(genome, rng, 20)
    write_fastq(os.path.join(reads, "enz"), pairs)
    with open(os.path.join(reads, "enz_truth.tsv"), "w") as fh:
        for row in truth:
            fh.write("\t".join(map(str, row)) + "\n")
    for label, bed in (("A", bed_a), ("B", bed_b)):
        with open(os.path.join(reads, f"enz_{label}.bed"), "w") as fh:
            for row in bed:
                fh.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
