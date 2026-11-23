# test/unit/allele/

End-to-end tests of the allele / Hi-C-matrix tools in `src/main/java/edu/mit/compbio/bisulfitehic/utils/`,
run on tiny synthetic BAM / VCF / BED files that the tests build themselves (no real data needed).

Run:

```
bash test/unit/allele/run.sh [scratch_dir]
```

The tool classes are compiled from the current `src/` tree (not taken from the jar in `target/`);
classes and data go to `scratch_dir` (default: a new `mktemp -d` directory), nothing is written into the repository.
Exit status 0 = all passed, 1 = a check failed, 2 = compilation failed.

| File | Purpose | How generated |
|---|---|---|
| `run.sh` | Compiles the tests plus the tool sources they use, runs every test class | Written by hand, 2026-11 |
| `Synth.java` | Shared builder: periodic synthetic genome (one CpG per 20 bp, A/T SNP sites), C->T-strand reads with MD tag and TLEN=0, indexed BAM, tabix-indexed phased VCF | Written by hand, 2026-11 |
| `LongRangeAsmTest.java` | `LongRangeAsm`: soft-clipped read assigned to the right allele; per-CpG files not repeated across chromosomes; distance filter from coordinates (TLEN=0 pairs); `-diffChr` (trans mate labelled with its own contig, short cis pair excluded); missing mate record; header labels | Written by hand, 2026-11 |
| `AlleImbalanceExprTest.java` | `AlleImbalanceExpr`: soft-clipped read; genes with identical coordinates; BED 0-based start; p-value not printed as 0.000; header/row column counts (per SNP and `-mergedInSameGene`) | Written by hand, 2026-11 |
| `SepAlleleMethySCTest.java` | `SepAlleleMethySC`: CpG of a soft-clipped read is not mistaken for the SNP position | Written by hand, 2026-11 |
| `SepAlleleMethySC2Test.java` | `SepAlleleMethySC2`: mates matched by read name (singleton in the BAM); SNP-position skip applied only to the read carrying the SNP; soft-clipped read | Written by hand, 2026-11 |
| `EnhanceHicMatrixTest.java` | `EnhanceHicMatrix -bedpe`: anchor as wide as `-window` is analysed; reads without a CpG of `cpg_list.bed` are not used | Written by hand, 2026-11 |

All synthetic reads are on the C->T bisulfite strand, so these tests do not cover the G->A-strand
coordinate handling of `BisulfiteHicUtils`.
