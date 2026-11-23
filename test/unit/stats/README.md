# test/unit/stats

Unit and small end-to-end tests for the statistics classes in
`src/main/java/edu/mit/compbio/bisulfitehic/utils/`: `FisherExactTest`, `PoissonBinomial2`,
`SingleCellTest`, `TestIndepPoiBinom`, `TestIndepPoiBinomOnBam`, `SingleBaseResHicMatrix`,
`SingleBaseResHicMatrixV2`.

Run:

```
bash test/unit/stats/run.sh [scratchDir]
```

The script compiles the main sources and the tests into `scratchDir` (default: a new `mktemp`
directory, never `target/`) and runs each test in its own JVM. Exit code 0 means every check
passed. All inputs (text matrices, SAM/BAM, VCF) are synthetic and written by the tests into
`scratchDir/work`; nothing here needs a reference genome. Set `MAIN_CLASSES=<dir>` to test an
already compiled build (for example a pre-fix build) instead of the source tree.

| File | Purpose | How generated |
|---|---|---|
| `run.sh` | Compiles and runs all tests in this folder | Written by hand, 2026-11 |
| `T.java` | Assertion and file helpers shared by the tests (no test framework) | Written by hand |
| `StatsSynth.java` | Builds synthetic bisulfite read pairs (reference carried by the MD tag), an indexed BAM, a SAM and a tabix-indexed het-SNP VCF | Written by hand |
| `FisherExactTestTest.java` | One- and two-sided Fisher p-values against exact integer arithmetic for every 2x2 table with n <= 30 (ties in the two-sided test) | Written by hand |
| `PoissonBinomial2Test.java` | Poisson-binomial PMF/CDF against a direct convolution, including odd trial counts and p = 1 | Written by hand |
| `SingleCellTestOutputTest.java` | Output format of `SingleCellTest` (numeric row label, aligned NA columns) and `-regions` | Written by hand |
| `TestIndepPoiBinomWindowTest.java` | Output shape of `TestIndepPoiBinom`, including the last window when `-aveWindow` > 1 | Written by hand |
| `TestIndepPoiBinomOnBamTest.java` | `TestIndepPoiBinomOnBam` on two synthetic cells: one slot per cell and allele, and the permutation null | Written by hand |
| `SingleBaseResHicMatrixTest.java` | `SingleBaseResHicMatrix` (`v1`) and `SingleBaseResHicMatrixV2` (`v2`): concordance counts, read with a deletion, last CpG, `-permutations 0`, `-bedpe -skipProcessFirstRow`, `-regions` | Written by hand |
| `FILES.md` | One-line description of every file in this folder | Written by hand |

Notes:

- The CpG list in `SingleBaseResHicMatrixTest` uses the 1-based position of the C as the BED
  start, because that is the position both tools look up. Whether the tools should accept a
  standard 0-based BED instead is an open question and is not tested here.
- `TestIndepPoiBinomOnBamTest` checks the permutation score with a loose bound (< 1.5) because
  the permutations are not seeded.
