# test/regression/

End-to-end regression tests for the Bhmem aligner and helper scripts. They build a
tiny synthetic bisulfite reference (two 60 kb random chromosomes), simulate
directional and PBAT Hi-C read pairs, run Bhmem from the built jar and check the
output BAM. Runtime is about one minute; no external data is needed.

## Run

```
bash test/regression/run_tests.sh
```

Requires `java`, `samtools`, `python3`, and the bundled `bwa` + `libbwajni.so`
that `install.sh` builds under `jbwa/jbwa-1.0.0/`. On Quest, `java` is loaded from
the module `java/jdk-17.0.2+8` if it is not on `PATH`. Environment overrides
(`JAR`, `WORK`, `JAVA_MODULE`, `PYSAM_PYTHON`) are documented at the top of the
script. Exit status is 0 only if every test passes.

## What is tested

| Test | Checks |
|---|---|
| 1. directional library | QUAL is reversed exactly when SEQ is reverse-complemented; aligned bases match the reference apart from bisulfite conversions; the strand flag agrees with SEQ |
| 2. PBAT library (`-pbat -nonDirectional`) | Same three checks with the PBAT strand convention (flag inverted, SEQ on the reference forward strand) |
| 3. `-enzymeList` | A split read whose two placements tie on mapping quality is placed at the locus supported by restriction sites, for two opposite site lists |
| 4. sequence dictionary | Reference under a directory whose name matches the regex `.fa`; reference with a `.fai` but no `.dict`; reference given as a bare file name |
| 5. `mh_hic_filter.py` | Long-range cis pairs are counted from mate positions even when TLEN is 0 |

Tests 1 and 2 include reads that span a ligation junction, so supplementary
alignments are exercised.

## Files

| File | Purpose | How generated |
|---|---|---|
| `run_tests.sh` | Test runner: builds the reference index, runs Bhmem, calls the checkers | Written by hand, 2026-11-21 |
| `make_test_data.py` | Seeded generator for the reference, bisulfite genomes, FASTQs, enzyme-site lists | Written by hand, 2026-11-21 |
| `check_bam.py` | SEQ / QUAL / strand consistency checker; reads `samtools view` text, accepts plain or gzipped FASTQ, so it also works on real data (`--pbat`, `--max-mismatch-frac`) | Written by hand, 2026-11-21 |
| `check_enzyme.py` | Checker for test 3 | Written by hand, 2026-11-21 |
| `work/` | Scratch output of the last run (git-ignored, recreated on each run) | `run_tests.sh` |
