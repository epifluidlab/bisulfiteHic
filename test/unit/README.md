# test/unit/

Unit tests for the analysis classes and helper scripts, added 2026-11-23 together with
the fixes they cover. Every test builds its own tiny synthetic input (SAM/BAM, VCF, text
matrices); no external data is needed. Each test was first run against the unfixed code
and failed there.

```
bash test/unit/run_all.sh
```

| Folder | Covers | Run alone |
|---|---|---|
| `concordance/` | `BisulfiteHicUtils`, `MethyCorAcrossHiccups`, `MethyCorGlobal`, `MethyCorEachChrScBw` | `bash test/unit/concordance/run.sh` |
| `allele/` | `LongRangeAsm`, `AlleImbalanceExpr`, `SepAlleleMethySC`, `SepAlleleMethySC2`, `EnhanceHicMatrix` | `bash test/unit/allele/run.sh` |
| `stats/` | `FisherExactTest`, `PoissonBinomial2`, `SingleCellTest`, `TestIndepPoiBinom`, `TestIndepPoiBinomOnBam`, `SingleBaseResHicMatrix` | `bash test/unit/stats/run.sh` |
| `perl/` | `src/perl/demultiplex_ecker_scWGBS.pl` | `bash test/unit/perl/run.sh` |
| `python/` | `src/python/*.py`, `src/python/methy_sc_concordance/*.py` | `python -B -m unittest discover -s test/unit/python` (needs numpy, scipy, pysam, bx-python) |

The Java groups compile the sources under `src/main/java` into a temporary directory and
put those classes ahead of the jar on the classpath, so they test the working tree, not
the built jar, and never write into `target/`. They need a JDK 9 or newer (`--release 8`);
`JAVA_HOME` defaults to the Quest OpenJDK 11. The dependencies are taken from
`target/bisulfitehic-0.38-jar-with-dependencies.jar`, `lib/` and `jbwa/`, so the project
must have been built once.

The aligner itself (Bhmem) is covered by the end-to-end tests in `test/regression/`.
