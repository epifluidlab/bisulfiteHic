# docs/

Usage documentation for the package. Start with `USAGE.md`.

| File | What it is | How it was generated |
|---|---|---|
| `USAGE.md` | Entry point: install and run, reference preparation, the standard workflow, how to read a Bhmem BAM, which tool answers which question, shared conventions, what to check before trusting a result | Written by hand, 2026-11-26 |
| `aligners.md` | Reference for `Bhmem`, `BhmemMultRef`, `BsmemMultRef`, the installer scripts and the `bisulfitehicMap` wrapper | Written from a read-through of the source at commit `289ba88`, 2026-11-26 |
| `concordance.md` | Reference for the five `MethyCor*` tools | As above |
| `allele_specific.md` | Reference for `LongRangeAsm`, `AlleImbalanceExpr`, `SepAlleleMethySC`, `SepAlleleMethySC2`, `EnhanceHicMatrix` | As above |
| `matrix_statistics.md` | Reference for `SingleBaseResHicMatrix`, `SingleBaseResHicMatrixV2`, `SingleCellTest`, `TestIndepPoiBinom`, `TestIndepPoiBinomOnBam` | As above |
| `scripts.md` | Reference for the Python, Perl and R helper scripts | As above; the scripts were also run on small synthetic inputs |

How the pages were checked:

- The option table of every Java tool was compared with the `@Option` declarations in the source: same options, none missing, none extra.
- Every Java tool was started without arguments with the documented command shape and printed its usage.
- For `Bhmem`, the pair-selection order, the `-pbat` flag handling, the TLEN reset and the effect of `samtools fixmate` were checked against the code and the regression BAMs.
- Output column lists were spot-checked against the code for `MethyCorAcrossHiccups`; the others were taken from the read-through and are not individually re-checked.

When a tool changes, update its section and its entry under "Known limitations" on the same page.
