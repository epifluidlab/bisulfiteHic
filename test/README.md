# test/

| File | Purpose | How generated |
|---|---|---|
| `test_R1.fq.gz`, `test_R2.fq.gz` | Small example read pair for trying the mapping commands in the top-level README (needs a real bisulfite-indexed reference) | Shipped with the original bisulfitehic repository |
| `regression/` | Self-contained end-to-end regression tests of the Bhmem aligner (synthetic reference and reads). See `regression/README.md` | Added 2026-11-21 |
| `unit/` | Unit tests of the analysis classes and helper scripts on synthetic inputs (`bash test/unit/run_all.sh`). See `unit/README.md` | Added 2026-11-23 |
