# test/unit/concordance

Unit tests for the methylation-concordance classes in
`src/main/java/edu/mit/compbio/bisulfitehic/utils/`:
`BisulfiteHicUtils`, `MethyCorAcrossHiccups`, `MethyCorGlobal`, `MethyCorEachChrScBw`.
(`MethyCorAcrossHiccupsSC` is not covered: it needs bigWig inputs.)

No real genomic data is used. Every BAM is built in memory by the test from a 31 bp synthetic
reference segment (reads carry an MD tag, so no reference FASTA is needed) and written to a
scratch directory.

## Files

| File | What it is | How it was made |
|---|---|---|
| `ConcordanceTest.java` | Plain Java class with a `main(scenario, workDir)`; exits 1 when a check fails. | Hand-written, 2026-11-21. |
| `run.sh` | Compiles the main sources and the test into a scratch directory (never `target/`) and runs each scenario in its own JVM. | Hand-written, 2026-11-21. |

## Running

```bash
bash test/unit/concordance/run.sh [scratchDir]
```

- `JAVA_HOME` selects the JDK (default: the Quest OpenJDK 11 install).
- `MAIN_CLASSES=<dir>` tests already compiled classes instead of compiling the source tree,
  for example a build of an older revision.

## Scenarios

| Scenario | What it checks |
|---|---|
| `helpers` | Pattern matcher ignores non-ACGT reference bytes (inserted bases, N); trailing insertion keeps reference and read the same length; per-location offsets are reference-forward for G->A reads. |
| `random` | `MethyCorAcrossHiccups -randomHiccupsInput` writes rows with the default `-randomSampleReads`. |
| `count` | `-outputCount` total-read-pair column is the same on every row. |
| `diffchr` | `-diffChr` rows are not dropped by the same-chromosome distance filter. |
| `short` | `-onlyShort` / `-shortLowerBound` use the mate distance, not TLEN. |
| `selfpair` | In `-onlyShort` mode a read overlapping both anchors is not paired with itself. |
| `unpaired` | `-useUnpaired` reads reach the methylation caller. |
| `global` | `MethyCorGlobal -useBadMate` uses non-proper pairs. |
| `oddbws` | `MethyCorEachChrScBw -shuffleWithinCell` rejects an odd number of `-bws`. |
