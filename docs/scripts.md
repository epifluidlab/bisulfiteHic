# Helper scripts (Python, Perl, R)

Part of the [usage guide](USAGE.md). Paths are relative to the repository root.

- Python packages imported across the scripts: `pysam`, `numpy`, `scipy`, `bx-python`, `pyBigWig`, `tqdm`, and `fastbetabino` (only `betabinom_test_mgch.py`). The R script needs `tidyverse`, `TopDom`, `strawr`, `hicrep` and `logging`.
- Statements marked "(run check)" were confirmed by running the script on small synthetic inputs. The Python unit tests (`test/unit/python/`, 28 tests) and the Perl test pass on Python 3.12, pysam 0.22.
- `eigen_adjust_direction.py` and `call_topdom_on_hic.R` were not run; their sections come from the code only.
- Sentences that start with "Not verified:" were not checked.

Contents: [Hi-C conversion and QC](#hi-c-conversion-and-qc-scripts) · [GCH/HCG methylation tracks](#gchhcg-methylation-track-scripts) · [Single-cell concordance](#single-cell-methylation-concordance-scripts-srcpythonmethy_sc_concordance) · [Perl and R](#perl-and-r-scripts) · [Known limitations](#known-limitations)

---

Which script is current, in short:

| Task | Current | Legacy |
|---|---|---|
| BAM to Juicer text | `src/python/sam2juicer_new.py` (Python 3; used by README Analysis 3 and by the sc_NOMeHiC pipeline) | `src/python/sam2juicer.py` (Python 2 only) |
| Read/pair QC summary | `src/python/mh_reads_summary.v2.py` (per fragment; the one the pipeline calls) | `src/python/mh_reads_summary.py` (per read; still maintained and the only one with a unit test) |

Invocations in the sc_NOMeHiC pipeline (Snakemake rules):

```
# rule demultiplex_fastqc_trim
perl {params.bisulfitehic}/src/perl/demultiplex_ecker_scWGBS.pl \
{params.out_prefix} {input.index_file} {input.r1} {input.r2} \
2> {log}

# rule qc
python {params.bisulfitehic}/src/python/mh_reads_summary.v2.py \
--in_cram {input.sorted_bam} --out_summary {params.summary_tmp} \
2> {log}
gzip -f {params.summary_tmp}

# rule scbam2hic
samtools view -@ {threads} -bh -q 30 -f 1 -F 1804 {input.calmd_bam} > {output.good_reads} && \
python {params.bisulfitehic}/src/python/sam2juicer_new.py \
-s {output.good_reads} -f {params.restriction_sites} 2> {log} | gzip -nc > {output.hic}
```

In those rules `params.out_prefix` is `01.demul_fastq_snakemake/{prefix}`, `input.r1`/`input.r2` are
`{prefix}_R1_001.fastq.gz`/`{prefix}_R2_001.fastq.gz`, the demultiplexed outputs are
`01.demul_fastq_snakemake/{prefix}.{index}.R1.fastq.gz` and `.R2.fastq.gz`, and `input.sorted_bam` is
`04.alignment_snakemake/{prefix}.{index}_sorted_by_name.calmd.bam`, which is produced with `samtools sort -n`.

---

## Hi-C conversion and QC scripts

### `src/python/sam2juicer_new.py`

**Purpose** — Converts paired alignments in a SAM or BAM file into one text line per read pair in the 16-field Juicer text format, written to standard output. This is the current converter; README Analysis 3 and the pipeline use it.

**Command** — `python sam2juicer_new.py -s <sam_or_bam> -f <restriction_site_file>`. It needs Python 3.6 or newer (it uses an f-string) and the `pysam` package; the other imports (`sys`, `argparse`, `bisect`) are standard library.

| Option | Default | Meaning |
|---|---|---|
| `-s`, `--sam` | required | Input alignment file. |
| `-f`, `--frag` | required | Restriction-site ("fragment site") file. |

**Inputs**
- Alignment file. It is opened with `pysam.AlignmentFile(path, "r")`. A BAM file is read correctly with this mode (run check; the pipeline and the README both pass a BAM). No index is needed. No sort order is needed: the two records of a pair are matched by read name through an in-memory dictionary, in whatever order they appear.
- Restriction-site file. Plain text, one line per contig, fields separated by a single space: the contig name followed by the site positions as integers, for example `chr1 1000 5000 100000`. There is no header. Contig names are normalised when the file is loaded: every `chr` is removed and `MT` is replaced by `M`.

**Output** — Standard output, one line per pair, 16 fields separated by single spaces, no header:

| # | Field | How it is filled |
|---|---|---|
| 1 | strand1 | `flag & 16` of the first end: `0` for forward, `16` for reverse. |
| 2 | chr1 | Reference name of the first end exactly as in the BAM (not normalised). |
| 3 | pos1 | `reference_start` of the first end: the 0-based leftmost aligned position, also for reverse-strand reads. |
| 4 | frag1 | `bisect(sites, pos1)`: the number of listed sites on that contig that are less than or equal to pos1. It is 0 for a read left of the first site. |
| 5 | strand2 | As field 1, for the second end. |
| 6 | chr2 | As field 2, for the second end. |
| 7 | pos2 | As field 3, for the second end. |
| 8 | frag2 | As field 4, for the second end. |
| 9 | mapq1 | `mapping_quality` of the first end. |
| 10 | cigar1 | CIGAR string of the first end. |
| 11 | sequence1 | `query_sequence` of the first end as stored in the file. |
| 12 | mapq2 | `mapping_quality` of the second end. |
| 13 | cigar2 | CIGAR string of the second end. |
| 14 | sequence2 | `query_sequence` of the second end. |
| 15 | readname1 | Read name of the first end. |
| 16 | readname2 | Read name of the second end. |

Which record is "first": the record on the contig with the smaller index in the BAM header; on the same contig, the one with the smaller fragment index; then the one with the smaller strand value (forward before reverse); then the one with the smaller position. If everything ties, the record that appears later in the file is written first.

Example (run check): a pair with one end at SAM `chr1:501` forward, MAPQ 60, and the other at `chr2:3001` reverse, MAPQ 37, gives `0 chr1 500 0 16 chr2 3000 1 60 10M CCCCCCCCCC 37 10M AAAAAAAAAA r1 r1`.

Not verified: this field order is the 16-field "long format" accepted by `juicer_tools pre` (this comes from the Juicer documentation, not from this repository; the repository only shows that the README feeds the output to `juicer_tools pre -d -q 30 -f <sites>`).

**Behaviour notes**
- The script applies no filter at all: no MAPQ filter, no duplicate filter, no proper-pair check, no secondary or supplementary check. Any filtering has to be done before (the pipeline uses `samtools view -q 30 -f 1 -F 1804`) or after (README: `juicer_tools pre -q 30`).
- For the fragment lookup the BAM contig name is normalised the same way as the site-file names (`chr` removed, `MT` to `M`), so `MT`, `chrM` and `M` all match each other.
- If a contig of a pair is not in the site file the script stops with exit status 1 and a one-line message on stderr: `Contig '<name>' is missing from the fragment file ...`.
- Output order is the order in which pairs are completed in the input. The script does not sort.

**Pitfalls**
- An unmapped record without a reference name makes the script fail with `AttributeError: 'NoneType' object has no attribute 'replace'` (run check). Remove unmapped reads first.
- Only the first two records with the same read name are paired. A third record (for example a supplementary alignment, which `-F 1804` does not remove) is kept in memory and never written unless a fourth record with the same name appears (run check).
- Reads whose mate never appears stay in memory until the end, so memory use grows with the number of unpaired records.
- Positions are 0-based (SAM POS minus 1).
- The site file is split on single spaces. Tabs or double spaces are not handled.

**Minimal example** (README Analysis 3)
```
python src/python/sam2juicer_new.py -s sample.good_reads.bam -f /path/to/restriction_sites/hg19_DpnII.txt | gzip -fc > sample.hic.txt.gz
```

### `src/python/sam2juicer.py`

**Purpose** — The original Python 2 version of the converter (header comment: Yunjiang Qiu, 2016). It produces the same 16 fields as `sam2juicer_new.py`. It is not the current script.

**Command** — `python2 sam2juicer.py -s <sam_or_bam> -f <restriction_site_file>`. It uses the Python 2 iterator call `bam.next()` and the old pysam names (`pysam.Samfile`, `getrname`, `qname`, `pos`, `mapq`, `seq`, `tid`). Under Python 3.12 with pysam 0.22.1 it stops immediately with `AttributeError: 'pysam.libcsamfile.Samfile' object has no attribute 'next'` (run check). Not verified: that it still runs correctly under Python 2 with an old pysam (not tested here).

| Option | Default | Meaning |
|---|---|---|
| `-s`, `--sam` | required | Input alignment file. |
| `-f`, `--frag` | required | Restriction-site file. |

**Inputs** — Same as `sam2juicer_new.py`.

**Output** — Same 16 space-separated fields, in the same order and with the same ordering rule, on standard output.

**Behaviour notes** — How it differs from `sam2juicer_new.py`:
1. It needs Python 2 and the old pysam API; the new script uses the Python 3 iterator and current pysam attribute names.
2. For the fragment lookup it removes `chr` from the BAM contig name but does not map `MT` to `M`, while the site-file keys do have `MT` mapped to `M`. A BAM contig named `MT` is therefore never found. The new script normalises both sides the same way.
3. Any error in the fragment lookup is caught by a bare `except` and reported as `Something wrong locating fragments for <a> or <b>.`. The new script catches only a missing contig and names it.
4. The read loop treats any exception from `bam.next()` as end of file, so a read error in the middle of a file ends the output without an error message.
5. The first record is fetched before the loop, so an input without records raises an uncaught `StopIteration`.

**Pitfalls** — Do not use it in a Python 3 environment. The points listed for `sam2juicer_new.py` (no filtering, 0-based positions, unmapped reads, third records) apply here as well.

**Minimal example**
```
python2 src/python/sam2juicer.py -s sample.bam -f /path/to/restriction_sites/hg19_DpnII.txt > sample.hic.txt
```

### `src/python/mh_hic_filter.py`

**Purpose** — Counts high-quality, non-duplicate, long-range cis read pairs in a name-sorted alignment file and sets a 0/1 flag that tells whether any chromosome has too few long-range cis contacts for its length. It is meant as a per-cell quality filter.

**Command** — `python mh_hic_filter.py --in_cram <name_sorted.bam> --out_summary <out.txt>`. Python 3 with `pysam` (run check on Python 3.12.8, pysam 0.22.1).

| Option | Default | Meaning |
|---|---|---|
| `--in_cram` | none (must be given) | Input alignment file, sorted by read name. |
| `--out_summary` | none (must be given) | Output text file. |

Neither option is declared `required` in argparse; leaving one out gives a Python error instead of a usage message.

**Inputs** — The file is opened with pysam mode `"rc"` (CRAM). A BAM file is read correctly as well (run check). Not verified: real CRAM input is not exercised by any test. The two ends of a pair must be adjacent records with identical read names (no `/1`, `/2` suffix), which means name-sorted. No index is needed. The header `@SQ` lines (`SN`, `LN`) are used.

**Output** — A text file with four lines, each `label<TAB>integer`:

| Line | Rule |
|---|---|
| `UniqMappedMapQ30NoPcrCisMore1kb:` | Kept pairs (see below) whose two starts are at least 1,000 bp apart. |
| `UniqMappedMapQ30NoPcrCisMore10kb:` | The same with at least 10,000 bp. |
| `UniqMappedMapQ30NoPcrCisMore20kb:` | The same with at least 20,000 bp. |
| `PerChromCisLongFilter:` | `1` if at least one considered contig has a number of 10 kb-or-longer cis pairs that is less than or equal to its length in Mb (rounded down); otherwise `0`. |

**Behaviour notes**
- A pair is two consecutive records with the same read name. It is kept when the second record is flagged paired, both records are mapped, neither is secondary, both have MAPQ of at least 30, the second record is not flagged duplicate, and the second record's reference equals its mate reference (cis).
- The proper-pair flag is not consulted. Supplementary and QC-fail flags are not consulted.
- The cis distance is `abs(reference_start of record 2 - reference_start of record 1)`, that is, from the two leftmost aligned positions. TLEN is not used.
- Contigs considered for the per-chromosome flag are hard-coded: every `@SQ` contig except names that start with `chrUn`, end with `random`, or equal `chrM`, `chrX`, `chrY`.

**Pitfalls**
- The exclusion list only recognises UCSC-style names. Any other extra contig (a spike-in such as `lambda`, `chrEBV`, alt or decoy contigs, or `X`/`Y`/`MT` in a genome without the `chr` prefix) is included in the per-chromosome test. A contig shorter than 1 Mb has threshold 0, so it sets the flag to 1 whenever it has no 10 kb cis pair. With a header of `chr1`, `chrX`, `lambda` and one long cis pair on `chr1`, the flag is 1; without the `lambda` line it is 0 (run check).
- The remembered first record is not replaced after a match. If three records share a read name, the third is compared with the first again and is counted as one more pair (run check: a pair plus one supplementary record gives counts of 2).
- Only the second record's duplicate flag and paired flag are tested.
- Input that is not name-sorted gives counts near zero without any warning, because mates are then not adjacent.

**Minimal example** (built from the argparse definition; there is no test, README or pipeline example for this script)
```
python src/python/mh_hic_filter.py --in_cram cell_sorted_by_name.calmd.bam --out_summary cell.hic_filter.txt
```

### `src/python/mh_reads_summary.py`

**Purpose** — Per-read mapping and Hi-C contact statistics (version 1). Every alignment record is counted on its own, so a pair contributes two counts. It is not the script the pipeline calls, but it was updated on 2026-11-23 and it is the only reads-summary script with a unit test.

**Command** — `python mh_reads_summary.py --in_cram <bam> --out_summary <out.txt>`. Python 3 with `pysam`.

| Option | Default | Meaning |
|---|---|---|
| `--in_cram` | none (must be given) | Input alignment file. |
| `--out_summary` | none (must be given) | Output text file. |

**Inputs** — Opened with pysam mode `"rc"`; a BAM file works (unit test and run check). No sort order and no index are needed, because each record is judged from its own flags and its own mate fields (RNEXT, PNEXT). The spike-in contig must be named exactly `lambda` to be recognised.

**Output** — A text file with 18 lines, each `label<TAB>count<TAB>percent%`. The percentage is an unrounded Python float (for example `66.66666666666667%`); it is `0.0%` when the denominator is 0. The first line always shows the literal `100%`. Each category is a subset of the one it is nested under:

| # | Label | Rule (counts records) | Percent of |
|---|---|---|---|
| 1 | `TotalReads:` | Every record in the file. | literal `100%` |
| 2 | `TotalMapped:` | Not unmapped. | TotalReads |
| 3 | `UniqMapped:` | Mapped and not secondary. | TotalReads |
| 4 | `UniqMappedToLambda:` | UniqMapped on contig `lambda`. | TotalReads |
| 5 | `UniqPairedToLambda:` | Line 4, and paired, mate mapped, mate contig `lambda`. | TotalReads |
| 6 | `UniqMappedToTargetSpecies:` | UniqMapped on any contig other than `lambda`. | TotalReads |
| 7 | `UniqPairedToTargetSpecies:` | Line 6, and paired with a mapped mate. | TotalReads |
| 8 | `UniqMappedNoPcr:` | Line 7, and not flagged duplicate. | line 7 |
| 9 | `UniqMappedNoPcrMinMapQ1:` | Line 8, and MAPQ of this record at least 1. | line 7 |
| 10 | `UniqMappedNoPcrMinMapQ1Trans:` | Line 9, and mate on a different contig. | line 9 |
| 11 | `UniqMappedNoPcrMinMapQ1Cis:` | Line 9, and mate on the same contig. | line 9 |
| 12 | `UniqMappedNoPcrMinMapQ1CisMore1kb:` | Line 11, and distance at least 1,000 bp. | line 11 |
| 13 | `UniqMappedNoPcrMinMapQ1CisMore20kb:` | Line 11, and distance at least 20,000 bp. | line 11 |
| 14 | `UniqMappedNoPcrMinMapQ30:` | Line 8, and MAPQ of this record at least 30. | line 7 |
| 15 | `UniqMappedNoPcrMinMapQ30Trans:` | Line 14, and mate on a different contig. | line 14 |
| 16 | `UniqMappedNoPcrMinMapQ30Cis:` | Line 14, and mate on the same contig. | line 14 |
| 17 | `UniqMappedNoPcrMinMapQ30CisMore1kb:` | Line 16, and distance at least 1,000 bp. | line 16 |
| 18 | `UniqMappedNoPcrMinMapQ30CisMore20kb:` | Line 16, and distance at least 20,000 bp. | line 16 |

**Behaviour notes**
- Cis distance is `abs(reference_start - next_reference_start)`, taken from the record's own position and its mate-position field. TLEN is not used.
- "Uniq" means "not secondary". Supplementary records are counted like primary ones.
- Only the MAPQ of the record itself is tested; the mate's MAPQ is not.
- The proper-pair flag is not consulted.

**Pitfalls**
- All counts are reads, not pairs. A long-range pair adds 2 to the cis lines (unit test).
- The numbers are not comparable with those of `mh_reads_summary.v2.py` (different unit, different order of filters, different denominators).

**Minimal example** (same arguments as the function call in `test/unit/python/test_hic_scripts.py`)
```
python src/python/mh_reads_summary.py --in_cram sample.bam --out_summary sample.summary.txt
```

### `src/python/mh_reads_summary.v2.py`

**Purpose** — Per-fragment (read pair) mapping and Hi-C contact statistics (version 2). This is the version the sc_NOMeHiC pipeline calls for per-cell QC.

**Command** — `python mh_reads_summary.v2.py --in_cram <name_sorted.bam> --out_summary <out.txt>`. Python 3 with `pysam`.

| Option | Default | Meaning |
|---|---|---|
| `--in_cram` | none (must be given) | Input alignment file, sorted by read name. |
| `--out_summary` | none (must be given) | Output text file (plain text; the pipeline gzips it afterwards). |

**Inputs** — Opened with pysam mode `"rc"`; a BAM works (the pipeline passes a BAM; run check). The two ends of a pair must be adjacent records with identical read names, which means name-sorted (`samtools sort -n`). No index is needed. The spike-in contig must be named exactly `lambda`.

**Output** — A text file with 13 lines, each `label<TAB>count<TAB>percent%` (unrounded float; `0.0%` when the denominator is 0; the first line shows the literal `100%`). "Record 2" is the second of two consecutive records with the same name, "record 1" the first:

| # | Label | Rule (counts pairs) | Percent of |
|---|---|---|---|
| 1 | `TotalFragments:` | Number of times a record has the same name as the remembered previous record. Read names that occur only once are not counted. | literal `100%` |
| 2 | `UniqMapped:` | Record 2 is flagged paired, both records are mapped, neither is secondary. | line 1 |
| 3 | `UniqMappedMapQ30:` | Line 2, and both records have MAPQ of at least 30. | line 1 |
| 4 | `UniqMappedMapQ30ToLambda:` | Line 3, and record 2 is on contig `lambda`. | line 3 |
| 5 | `UniqMappedMapQ30ToTargetSpecies:` | Line 3, and record 2 is not on `lambda`. | line 3 |
| 6 | `UniqMappedMapQ30NoPcr:` | Line 3, and record 2 is not flagged duplicate. | line 3 |
| 7 | `UniqMappedMapQ30NoPcrToLambda:` | Line 6, and record 2 is on `lambda`. | line 3 |
| 8 | `UniqMappedMapQ30NoPcrToTargetSpecies:` | Line 6, and record 2 is not on `lambda`. | line 3 |
| 9 | `UniqMappedMapQ30NoPcrTrans:` | Line 8, and record 2's reference differs from its mate reference. | line 8 |
| 10 | `UniqMappedMapQ30NoPcrCis:` | Line 8, and record 2's reference equals its mate reference. | line 8 |
| 11 | `UniqMappedMapQ30NoPcrCisMore1kb:` | Line 10, and distance at least 1,000 bp. | line 8 |
| 12 | `UniqMappedMapQ30NoPcrCisMore10kb:` | Line 10, and distance at least 10,000 bp. | line 8 |
| 13 | `UniqMappedMapQ30NoPcrCisMore20kb:` | Line 10, and distance at least 20,000 bp. | line 8 |

**Behaviour notes**
- Cis distance is `abs(reference_start of record 2 - reference_start of record 1)`. TLEN is not used.
- Lambda/target, duplicate, and cis/trans are all decided from record 2 only.
- The proper-pair flag, the supplementary flag and the QC-fail flag are not consulted.
- How it differs from `mh_reads_summary.py`: it counts pairs instead of reads; it needs name-sorted input; it requires MAPQ of at least 30 on both ends before any other category and has no MAPQ 1 tier; it splits lambda from target after the MAPQ filter; it adds a 10 kb category; the cis-distance lines are percentages of the non-duplicate target pairs, not of the cis pairs; the distance comes from the partner record instead of the mate-position field.

**Pitfalls**
- Input that is not name-sorted gives `TotalFragments` near zero without any warning.
- The remembered first record is not replaced after a match. A third record with the same name (for example a supplementary alignment) is compared with the first one again and counts as one more fragment in every line it passes (run check: a pair plus one supplementary record gives `TotalFragments: 2`).
- Reads without a mate record in the file are invisible: they are not in `TotalFragments`, so the percentages are relative to pairs present in the file, not to sequenced fragments.
- This file has no unit test; the test suite covers only `mh_reads_summary.py`.

**Minimal example** (pipeline Snakefile, rule `qc`)
```
python src/python/mh_reads_summary.v2.py --in_cram cell_sorted_by_name.calmd.bam --out_summary cell.summary.txt
```

### `src/python/interact_density_around_regions.py`

**Purpose** — A virtual-4C style profile: for reads that overlap one "from" region, it counts where their mates start, in fixed-size bins around the centre of a second region.

**Command** — `python interact_density_around_regions.py --in_bam <sorted_indexed.bam> --out_summary <out.txt> [options]`. Python 3 with `pysam` and `bx-python` (`bx.intervals.intersection`).

| Option | Default | Meaning |
|---|---|---|
| `--in_bam` | none (must be given) | Coordinate-sorted, indexed BAM. |
| `--out_summary` | none (must be given) | Output text file. |
| `--from_region` | `chr17:85611135-85828994` | Viewpoint region `chrom:start-end`; reads are fetched from it. |
| `--plot_center` | `chr1:115724771-115798945` | Region `chrom:start-end`; only its midpoint `int((start+end)/2)` is used. |
| `--plot_bin_size` | `50000` | Bin size in bp. |
| `--plot_boundary` | `5000000` | Half-width of the plotted window in bp. |
| `--min_mapq` | `30` | Minimum MAPQ of the fetched read. |
| `--min_insert_size` | `20000` | Minimum distance between read start and mate start for same-chromosome pairs. |

**Inputs** — The BAM is opened in mode `"rb"` and read with `fetch(contig, start, stop)`, so it must be coordinate-sorted and indexed. Regions are parsed by splitting on `:` and `-`; start and end are passed to pysam as given (0-based, end-exclusive). Mate position comes from the RNEXT and PNEXT fields of each fetched read.

**Output**
- Standard output: one line with the total number of counted reads.
- `--out_summary`: two lines, every value followed by a tab (so each line ends with a tab).
  - Line 1: bin labels, integers from `-(num_bins/2)` to `num_bins/2 - 1`, where `num_bins = int(2 * plot_boundary / plot_bin_size)` (200 labels, -100 to 99, with the defaults).
  - Line 2: the count for each bin. The bin with label `i` holds mates whose start is in `[centre + i*bin_size, centre + (i+1)*bin_size)`.

**Behaviour notes**
- A read is counted when it is paired, mapped, has a mapped mate, is not secondary, not duplicate, not QC-fail, not supplementary, has MAPQ of at least `--min_mapq`, and is either on a different chromosome than its mate or at least `--min_insert_size` away from it (`abs(reference_start - next_reference_start)`; TLEN is not used). The mate's MAPQ is not tested and the proper-pair flag is not consulted.
- The mate must be on the `--plot_center` chromosome and its start must fall in one of the bins. Bin index uses floor division, so the bin left of the centre is `[-bin_size, 0)`.
- The unit is reads: if both ends of a pair overlap the from-region and both mates fall in the window, the pair adds 2.

**Pitfalls**
- The default regions are specific loci left in the code; always pass both regions.
- Contig names that contain `:` or `-` cannot be parsed.
- If the bin size does not divide the boundary, mates in the partial bins at the edges are not counted (unit test: bin 30, boundary 100 gives 6 bins that cover -90 to +90).

**Minimal example** (values from `test/unit/python/test_hic_scripts.py`)
```
python src/python/interact_density_around_regions.py --in_bam sample.sorted.bam --out_summary density.txt --from_region chr2:1000-5000 --plot_center chr1:10000000-10000000 --plot_bin_size 50000 --plot_boundary 5000000
```

### `src/python/eigen_adjust_direction.py`

**Purpose** — Takes a per-bin eigenvector for one chromosome (one value per line, as written by the Juicer eigenvector tool), drops bins with low mappability, flips the sign if needed so that positive values are the gene-dense side, and writes a four-column bed-like file.

**Command** — `python eigen_adjust_direction.py --in_eig_file <eig.txt> --out_file <out.bed> --chrom <name> [options]`. Imports `numpy`, `pyBigWig` and `bx-python`. Not verified: this script was not executed in this review and has no unit test; everything below is from reading the code.

| Option | Default | Meaning |
|---|---|---|
| `--in_eig_file` | none (must be given) | Eigenvector file for one chromosome. |
| `--out_file` | none (must be given) | Output file. |
| `--chrom` | none (must be given) | Chromosome to process. |
| `--chrom_size` | `/home/dnaase/startup/data/genomes/hg19/human_g1k_v37.common_chrom.size` | Chromosome sizes file. |
| `--w_len` | `2500000` | Bin size in bp; it must equal the resolution of the eigenvector. |
| `--gene_density` | `/home/dnaase/startup/data/annotations/gene/gencode.v30.b37.gene_density.bedgraph` | Gene interval file. |
| `--mappability` | `/home/dnaase/startup/data/genomes/hg19/human_g1k_v37.51mer.mappability.bw` | Mappability bigWig. |
| `--min_mappability` | `0.75` | A bin is kept only if its mean mappability is greater than this value. |

**Inputs**
- `--in_eig_file`: one number per line, no header; line `k` (starting at 0) is the bin that starts at `k * w_len`.
- `--chrom_size`: tab-separated, exactly two columns per line (name, size); no header.
- `--gene_density`: tab-separated, exactly four columns per line (chrom, start, end, value); no header. Only the intervals are used: the density of a bin is the number of intervals that overlap it, times 1,000,000, divided by `w_len`. The fourth column is ignored.
- `--mappability`: bigWig; the mean over each bin is read with `pyBigWig.stats(..., exact=True)`.
- `--chrom` must be spelled the same way in all three files.

**Output**
- Standard output: two lines, `pos sum:<x>` and `neg sum:<y>`. Despite the labels these are the mean gene densities of the kept bins with eigenvalue at least 0 and below 0.
- `--out_file`: tab-separated, no header, one line per kept bin: chrom, bin start, bin end (`min(chromosome length, start + w_len)`), adjusted eigenvalue.

**Behaviour notes** — The sign of all values is flipped when the mean gene density of the non-negative bins is lower than that of the negative bins. Bins with mean mappability at or below the threshold, or without mappability data, are left out of the output.

**Pitfalls**
- The three default paths point to `/home/dnaase/...` (hg19/b37 files) and will not exist elsewhere; pass `--chrom_size`, `--gene_density` and `--mappability` explicitly.
- If the eigenvector file has fewer lines than the chromosome has bins the script fails with an index error; if it has more, the extra lines are silently ignored. Nothing checks that `--w_len` matches the eigenvector resolution.
- A `.fai` file cannot be used as `--chrom_size` because lines must have exactly two columns.
- Bins whose eigenvalue is NaN take no part in the sign decision but are still written.

**Minimal example** (built from the argparse definition; there is no test, README or pipeline example)
```
python src/python/eigen_adjust_direction.py --in_eig_file chr1.eigen_2.5Mb.txt --out_file chr1.eigen_adjusted.bed --chrom 1 --chrom_size genome.chrom.size --w_len 2500000 --gene_density gene_density.bedgraph --mappability mappability.bw
```

---

## GCH/HCG methylation track scripts

All three read the Bis-SNP "6plus2" bed: tab-separated `chrom, start, end, name, score, strand, methylation (percent, 0-100), coverage (reads)`. Column numbers below are 1-based.

### `src/python/smooth_mgch.py`

**Purpose** — Turns a per-site methylation track into a per-base track: between each site and the next one it fits a straight-line smoothing spline through the neighbouring sites and writes one value for every base.

**Command** — `python smooth_mgch.py --in_bed <in.bed.gz> --out_bed <out.bed[.gz]> --chrom_size <sizes> [options]`. Python 3 with `pysam`, `numpy`, `scipy` (`UnivariateSpline`) and `bx-python`.

| Option | Default | Meaning |
|---|---|---|
| `--in_bed` | none (must be given) | bgzip-compressed, tabix-indexed bed. |
| `--out_bed` | none (must be given) | Output; written with gzip if the name ends in `.gz`, else plain text. |
| `--chrom_size` | none (must be given) | Tab-separated file; column 1 name, column 2 size. |
| `--window` | `1000` | Step in bp by which the neighbourhood grows on each side. |
| `--max_window` | `10000` | Largest reach of the neighbourhood on each side. |
| `--min_gch` | `100` | The neighbourhood stops growing once it holds at least this many sites. |
| `--column` | `7` | 1-based column that holds the value to smooth. |
| `--s` | `1` | Smoothing factor passed to `UnivariateSpline(k=1, s=...)`; `0` gives exact linear interpolation. |

**Inputs** — `--in_bed` must be bgzip + tabix (`pysam.TabixFile`); columns 1-3 and the `--column` column are used. Chromosomes are processed in the order of the `--chrom_size` file, and only those present in the tabix index.

**Output**
- Standard output: `chrom:<name>` and `with size of GCH:<number of sites>` for each chromosome.
- `--out_bed`: tab-separated, no header, four columns: chrom, start, end, value rounded to 4 decimals.
  - For a site that has a following site within 2,500 bp: one line for every base from the site's start up to the base before the next site's start (`end = start + 1`); the value is the spline evaluated at that base, with negative values set to 0.
  - For a site without a following site within 2,500 bp (including the last site of a chromosome): one line with the site's own start, end and unsmoothed value.

**Behaviour notes**
- The 2,500 bp limit is not an option: it is the default `max_dist` of `bx` `Intersecter.after()`, which the script calls without arguments (run check: sites at 1000, 1010, 6000, 6004 give per-base lines for 1000-1009 and 6000-6003 and single lines for 1010 and 6004).
- The neighbourhood starts as `[start - window, end)` plus `[next_start, next_start + window)` and grows outwards by `--window` on both sides until it holds `--min_gch` sites or reaches `--max_window`, or until both sides hit the contig ends.
- The spline degree is fixed at 1. Values are not capped at 100.

**Pitfalls**
- The output is per base and therefore much larger than the input.
- A `.gz` output is plain gzip, not bgzip; it has to be re-compressed with bgzip before it can be tabix-indexed.
- Gaps longer than 2,500 bp between sites are not filled.

**Minimal example** (option values from `test/unit/python/test_mgch_scripts.py`)
```
python src/python/smooth_mgch.py --in_bed sample.GCH.6plus2.bed.gz --out_bed sample.GCH.smooth.bed.gz --chrom_size genome.chrom.size --window 1000 --max_window 10000 --min_gch 100 --column 7 --s 1
```

### `src/python/binom_test_mgch.py`

**Purpose** — For every site in a 6plus2 bed, tests with a one-sided binomial test whether the site is more methylated than its local neighbourhood, and appends the p-value.

**Command** — `python binom_test_mgch.py --in_bed <in.bed> --out_bed <out.bed> [options]`. Python 3 is required (the background rate is computed with true division of two integers). Imports `scipy.stats` and `bx-python`. It uses `scipy.stats.binomtest` when it exists and falls back to `binom_test` on older scipy.

| Option | Default | Meaning |
|---|---|---|
| `--in_bed` | none (must be given) | Plain-text 6plus2 bed. |
| `--out_bed` | none (must be given) | Plain-text output. |
| `--window` | `100` | Flank size in bp on each side, and the step by which the flanks grow. |
| `--min_gch` | `3` | Keep growing while fewer than this many neighbouring sites were found. |
| `--min_point` | `10` | Keep growing while the neighbours have less than this total coverage. |
| `--max_window` | `1000` | Largest flank reach in bp. |
| `--global_acc` | `0.285` | Background probability used only when the neighbourhood has zero coverage. |

**Inputs** — Uncompressed text read with `open()`; exactly 8 tab-separated columns per line. Lines that start with `track` or `#` are skipped. Score (column 5) and coverage (column 8) must be integers. The whole file is loaded into memory and read twice. No sorting or index is needed.

**Output** — `--out_bed`: every data line of the input followed by a tab and the p-value formatted with `%.6g` (9 columns). Skipped `track`/`#` lines are not copied. Standard output: `build interval tree with size:<n>`.

**Behaviour notes**
- Methylated reads of a site: `round(methylation * coverage / 100)`.
- Background: all other sites in `[start - window, start)` and `[end, end + window)`; while there are fewer than `--min_gch` sites or less than `--min_point` total coverage, the next ring of `--window` bp is added on each side, up to `--max_window`. The background probability is methylated reads divided by coverage summed over those sites. The site itself is not part of its background.
- Test: `P(X >= methylated | n = coverage, p = background)`, alternative `greater`. A site with coverage 0 gets p = 1.

**Pitfalls**
- Input must be uncompressed and have exactly 8 columns; a gzip file or extra columns make it fail.
- `--min_gch` and `--min_point` only control how far the neighbourhood grows. If they are still not met at `--max_window`, the test is done with whatever was found.

**Minimal example** (default option values, as in the unit test)
```
python src/python/binom_test_mgch.py --in_bed sample.GCH.6plus2.bed --out_bed sample.GCH.binom.bed --window 100 --min_gch 3 --min_point 10 --max_window 1000
```

### `src/python/betabinom_test_mgch.py`

**Purpose** — For every site in a tabix-indexed 6plus2 bed, compares its methylation with two local backgrounds and the chromosome-wide background, and appends differences, binomial `-log10(p)` values, optionally a beta-binomial `-log10(p)`, and z-scores.

**Command** — `python betabinom_test_mgch.py --in_bed <in.bed.gz> --out_bed <out.bed[.gz]> --chrom_size <sizes> [options]`. Python 3 is required (true division of integer counts). Imports `scipy.stats`, `numpy`, `pysam`, `bx-python` and `from fastbetabino import *`. The `fastbetabino` import is unconditional, so the package must be importable even with `--binomial_only True`; it must provide `fit_alpha_beta(coverages, methylated_counts)`. Not verified: behaviour with the real `fastbetabino` package (it is not installed in the `scnomehic` env; the unit tests use a stand-in that returns alpha = beta = 1).

| Option | Default | Meaning |
|---|---|---|
| `--in_bed` | none (must be given) | bgzip + tabix 6plus2 bed. |
| `--out_bed` | none (must be given) | Output; gzip if the name ends in `.gz`, else plain text. |
| `--chrom_size` | none (must be given) | Tab-separated; column 1 name, column 2 size. Only listed chromosomes are processed. |
| `--window1` | `2500` | Flank size in bp of local background 1 (each side). |
| `--window2` | `5000` | Flank size in bp of local background 2 (each side). |
| `--min_gch` | `10` | Minimum number of sites for a local background to be used in its binomial test; also the minimum number of sites for a chromosome to be processed. |
| `--min_point` | `100` | Minimum total coverage for a local background to be used in its binomial test. |
| `--binomial_only` | `False` | Takes a value (`True/False`, `t/f`, `yes/no`, `y/n`, `1/0`). `True` skips the beta-binomial fit and its column. |
| `--two_sided` | `False` | Takes a value. `False` gives the one-sided test `less`; `True` gives `two-sided`. |

**Inputs** — `--in_bed`: bgzip + tabix, exactly 8 tab-separated columns, no header. Chromosomes are processed in `--chrom_size` order if present in the tabix index.

**Output** — `--out_bed`, tab-separated, no header: the 8 input columns followed by

| Column | Name in code | Content |
|---|---|---|
| 9 | `local_fold1` | Site methylation fraction minus the pooled fraction (methylated reads / coverage) of window 1; minus the chromosome-wide fraction if window 1 has zero coverage. |
| 10 | `local_fold2` | The same for window 2. |
| 11 | `global_fold` | Site fraction minus the chromosome-wide pooled fraction. |
| 12 | `binomial_p_w1` | `-log10` of the binomial p-value against the window 1 fraction; the chromosome-wide fraction is used instead if window 1 has fewer than `--min_gch` sites or less than `--min_point` coverage. |
| 13 | `binomial_p_w2` | The same for window 2. |
| 14 | `binomial_p_global` | `-log10` of the binomial p-value against the chromosome-wide fraction. |
| 15 | `binomial_p_final` | Maximum of columns 12-14 (not rounded). |
| 16 | `beta_binom_p_global` | Only when `--binomial_only False`: `-log10` of the beta-binomial lower-tail probability `cdf(methylated, coverage, alpha, beta)` with alpha and beta fitted on the chromosome. |
| 16 or 17 | `zscore_w1` | (Site fraction minus mean of the per-site fractions in window 1) divided by their standard deviation; chromosome-wide mean and standard deviation if window 1 has zero coverage. |
| 17 or 18 | `zscore_w2` | The same for window 2. |
| 18 or 19 | `zscore_global` | The same with the chromosome-wide mean and standard deviation. |

A line has 18 columns with `--binomial_only True` and 19 with `--binomial_only False` (unit test). All added values except column 15 are rounded to 4 decimals. Standard output gets progress lines (`chrom:<name>`, `with size of GCH:<n>`, and the fitted `global alpha`/`global beta`).

**Behaviour notes**
- Methylated reads: `round(methylation * coverage / 100)`.
- A local window is `[start - window, start)` plus `[end, end + window)`, clipped to the contig; the site itself is not in it. The chromosome-wide background does include the site.
- The default binomial alternative is `less`: large values in columns 12-15 mark sites that are less methylated than the background.
- A chromosome with fewer than `--min_gch` sites is skipped and produces no output lines.
- With more than 1,000,000 sites on a chromosome the beta-binomial fit uses a random sample of 1,000,000 sites drawn without a fixed seed.
- The beta-binomial column is always the lower tail; `--two_sided` does not change it.

**Pitfalls**
- Columns 9-11 are differences of fractions, not log2 fold changes (the comment next to `--out_bed` says log2 fold change).
- A p-value of 0 gives `inf`; a window in which all sites have the same fraction gives a zero standard deviation and an `inf` or `nan` z-score.
- A `.gz` output is plain gzip, not bgzip.
- Boolean options need an explicit value (`--binomial_only True`); they are not flags.

**Minimal example** (command line used in `test/unit/python/test_mgch_scripts.py`)
```
python src/python/betabinom_test_mgch.py --in_bed sample.GCH.6plus2.bed.gz --out_bed sample.GCH.betabinom.bed --chrom_size genome.chrom.size --min_gch 5 --min_point 10 --binomial_only True --two_sided False
```

---

## Single-cell methylation concordance scripts (`src/python/methy_sc_concordance/`)

### `src/python/methy_sc_concordance/call_methy_concord_sc.py`

**Purpose** — From per-cell window tables it builds, for every window, a Poisson-binomial distribution from the per-cell methylation probabilities, and then for every ordered pair of windows gives the probability of seeing at least the second window's observed methylated count under the first window's distribution.

**Command** — `python call_methy_concord_sc.py <input_dir> '<glob_pattern>' <output.gz>`. Three positional arguments read from `sys.argv`; there are no options and no usage message. Python 3 with `numpy`. The `PoiBin` class is defined inside this file; it does not import `poibin.py`.

| Argument | Meaning |
|---|---|
| 1 `input_dir` | Directory with the per-cell files. The script changes its working directory to it. |
| 2 `glob_pattern` | File pattern evaluated with `glob.glob` inside `input_dir`. Quote it so that the shell does not expand it. |
| 3 `output` | Output file, always written with gzip. |

**Inputs** — One tab-separated text file per cell, no header, at least 6 columns. Used columns (1-based): column 1 window name in the form `chrom:start-end`; column 3 coverage (integer; the line is used only if it is greater than 0); column 4 methylated count (integer); column 6 methylation probability between 0 and 1. Columns 2 and 5 are ignored.

**Output** — gzip text, tab-separated, no header, 8 columns: window 1 chrom, start, end; window 2 chrom, start, end; a name `chrom1:start1:end1:start2:end2`; the p-value `P(X >= k)`, where X follows the Poisson-binomial distribution of window 1's per-cell probabilities and k is the sum of column 4 over the cells of window 2. Standard output: `building distribution...` and `output...`.

**Behaviour notes**
- A window is used only if more than 10 files have a line for it with coverage above 0.
- All ordered pairs of used windows are tested, including a window with itself and windows on different chromosomes.
- A pair is written only when window 2's count is smaller than the number of cells of window 1; other pairs are left out silently.

**Pitfalls**
- Because of the change of directory, a relative output path is created inside `input_dir`.
- The number of output lines grows with the square of the number of windows.
- The name column does not contain the chromosome of window 2.
- Probabilities outside 0-1 raise a `ValueError`.

**Minimal example** (arguments as in `test/unit/python/test_sc_concordance.py`)
```
python src/python/methy_sc_concordance/call_methy_concord_sc.py /path/to/cells 'cell_*.txt' /path/to/out.concord.txt.gz
```

### `src/python/methy_sc_concordance/cal_concord_sc_at_regions.py`

**Purpose** — For each pair of regions in a bedpe file, takes the mean signal of both regions in every cell's bigWig and reports the Pearson and Spearman correlation across cells, optionally with permutation p-values.

**Command** — `python cal_concord_sc_at_regions.py --methy_file_list <list.txt> --bedpe <pairs.bedpe> --out_summary <out.txt> [options]`. Python 3.6 or newer (f-strings) with `numpy`, `scipy`, `pyBigWig` and `tqdm`.

| Option | Default | Meaning |
|---|---|---|
| `--methy_file_list` | none (must be given) | Text file with one bigWig path per line (one per cell). |
| `--bedpe` | none (must be given) | Region pairs. |
| `--out_summary` | none (must be given) | Output text file. |
| `--threads` | `1` | Number of worker processes (`multiprocessing.Pool`). |
| `--min_data` | `10` | A pair is reported only if more than this many cells have data at both regions. |
| `--permutation` | `0` | Number of permutations; `0` turns permutation p-values off. |

**Inputs**
- bedpe: tab-separated; the first six columns are chrom1, start1, end1, chrom2, start2, end2; further columns are ignored. Lines that start with `#` or have fewer than 6 columns are skipped. Lines whose coordinates are not integers are skipped with a message on stderr.
- bigWig files: the mean over each region is read with `pyBigWig.stats(chrom, start, end, exact=True)`. A cell is used for a pair only if both means exist.

**Output** — `--out_summary`, tab-separated, no header, one line per reported pair in input order:

| Column | Content |
|---|---|
| 1-6 | chrom1, start1, end1, chrom2, start2, end2. |
| 7 | Name `chrom1:start1:end1:chrom2:start2:end2`. |
| 8 | Pearson r (3 decimals; 5 decimals when `--permutation` > 0). |
| 9 | Pearson p-value (`%.6g`). |
| 10 | Spearman rho (3 decimals; 5 decimals when `--permutation` > 0). |
| 11 | Spearman p-value (`%.6g`). |
| 12 | Number of cells used. |
| 13 | Only with `--permutation` > 0: Pearson permutation p-value. |
| 14 | Only with `--permutation` > 0: Spearman permutation p-value. |

A progress bar goes to stderr.

**Behaviour notes**
- Permutation p-value: region 1's values are shuffled across cells N times; p = (number of permuted correlations greater than or equal to the observed one + 1) / (N + 1). It is one-sided (tests for positive correlation) and never 0.
- A region that reaches or passes the chromosome end in the bigWig is truncated to the chromosome end.

**Pitfalls**
- If exactly one of the two regions ends at or beyond the chromosome end, the other region's mean is never computed and the cell is dropped; this happens in every cell, so such a pair is never reported (run check: `chr1 100 200 chr1 9900 10000` on a 10,000 bp contig returns nothing, while both regions inside, or both at the end, are reported).
- The threshold is strict: with `--min_data 10` a pair needs at least 11 cells (run check).
- All bigWig files are opened again for every bedpe line, and the whole bedpe file is read into memory.
- Blank lines in `--methy_file_list` are not filtered out and are passed on as empty file names.
- Pairs with too few cells are left out of the output, so the output can have fewer lines than the bedpe.

**Minimal example** (built from the argparse definition; the unit tests call `process_line` directly)
```
python src/python/methy_sc_concordance/cal_concord_sc_at_regions.py --methy_file_list cells.GCH.bw.list --bedpe loops.bedpe --out_summary loops.GCH_concordance.txt --threads 4 --min_data 10 --permutation 0
```

### `src/python/methy_sc_concordance/test_indep_poibin.py`

**Purpose** — An analysis script (not a unit test, despite its name) that, for every pair of CpGs in a matrix of per-allele or per-cell methylation values, compares the observed standard deviation of the averaged pair with the sum of the two individual standard deviations.

**Command** — `python test_indep_poibin.py <input.txt.gz> <output.txt.gz>`. Two positional arguments read from `sys.argv`; no options. Python 3 with `numpy`.

**Inputs** — gzip text, tab-separated, no header, one CpG per line, every line with the same number of columns. Column 3 (1-based) is the CpG position (integer). Columns 7 onwards are the values, one per allele or cell; `NA` is read as missing. Columns 1, 2, 4, 5, 6 are not used. Example line from the unit test: `chr1<TAB>99<TAB>100<TAB>.<TAB>.<TAB>+<TAB>0<TAB>1<TAB>0<TAB>1<TAB>NA<TAB>1`.

**Output** — gzip text, tab-separated, no header and no row names: a square matrix with one row and one column per input CpG, in input order. Entry (i, j) is `nanstd((values_i + values_j) / 2) / (nanstd(values_i) + nanstd(values_j))`. It is `NA` when the two positions are more than 3,000,000 bp apart or when the observed standard deviation is not a number.

**Behaviour notes** — The 3,000,000 bp limit is hard-coded (`max_dist`). The diagonal is computed like any other entry.

**Pitfalls**
- The chromosome column is not read; distance is judged from column 3 alone, so the input should hold one chromosome.
- Time and output size grow with the square of the number of CpGs.
- The text `NA` is replaced anywhere in the line before splitting.

**Minimal example** (arguments as in `test/unit/python/test_sc_concordance.py`)
```
python src/python/methy_sc_concordance/test_indep_poibin.py cpg_by_allele.chr19.txt.gz cpg_pair_sd_ratio.chr19.txt.gz
```

### `src/python/methy_sc_concordance/poibin.py`

**Purpose** — A library module, not a command-line tool: class `PoiBin` implements the Poisson-binomial distribution (sum of independent Bernoulli trials with different success probabilities) with a fast-Fourier-transform method.

**Command** — None; running the file does nothing (`if __name__ == "__main__": pass`). Import it. It needs `numpy`. API: `PoiBin(probabilities)` with a one-dimensional list or array of values in 0-1; `.pmf(k)` gives `P(X = k)`; `.cdf(k)` gives `P(X <= k)`; `.pval(k)` gives `P(X >= k)`. `k` is an integer (Python or numpy) or a list of integers between 0 and the number of probabilities.

**Inputs** — The probability array given to the constructor. A `ValueError` is raised for input that is not one-dimensional or has values outside 0-1.

**Output** — Return values only. `pval` of a single integer returns `1` for `k = 0`, else `1 - cdf(k - 1)` clipped to 0-1; for a list it returns a numpy array of `1 - cdf(k) + pmf(k)` clipped to 0-1.

**Behaviour notes** — The full probability table is computed once in the constructor. Machine epsilon is added to every table entry.

**Pitfalls**
- Because of the added epsilon and FFT round-off, tail probabilities far below about 1e-16 are not resolved (the unit tests require agreement with the exact binomial to 10 decimal places).
- `call_methy_concord_sc.py` does not use this module; it carries its own copy of the class, which differs in one place (see Known limitations).

**Minimal example** (run check, prints 0.3703...)
```
python -c "import sys; sys.path.insert(0, 'src/python/methy_sc_concordance'); from poibin import PoiBin; print(PoiBin([0.2]*20).pval(5))"
```

---

## Perl and R scripts

### `src/perl/demultiplex_ecker_scWGBS.pl`

**Purpose** — Splits one pair of FASTQ files into one pair per cell, using a 6-base barcode at the start of read 1, and removes the barcode from read 1.

**Command** — `perl demultiplex_ecker_scWGBS.pl <outputPrefix> <index_file.txt> <R1.fastq[.gz|.bz2]> <R2.fastq[.gz|.bz2]>`. Four positional arguments. There are no options: the usage text says `[Options]` and `Getopt::Long` is loaded, but `GetOptions` is never called. With fewer than four arguments it prints the usage to stderr and exits with status 1. It calls the external programs `zcat` or `bzcat` (for compressed input) and `gzip`.

**Inputs**
- Index file: one barcode per line, nothing else (no header, no extra columns). The whole line is the key, so each barcode must be exactly 6 bases.
- Barcode position: the first 6 bases of the read 1 sequence line. Matching is exact (hash lookup, no mismatches allowed).
- FASTQ: four lines per record, the same number of records in both files and in the same order (read names are not compared). Compression is decided from the read 1 file name only: a name containing `.bz2` is read through `bzcat`, a name containing `.gz` through `zcat`, anything else as plain text. Read 2 is opened the same way as read 1.

**Output**
- For every barcode in the index file: `<outputPrefix>.<barcode>.R1.fastq.gz` and `<outputPrefix>.<barcode>.R2.fastq.gz`. They are written as plain `.fastq` and compressed with `gzip -f` at the end. A barcode that matches no read gives two empty gzip files.
- Read 1 records: header line unchanged, sequence without its first 6 bases, third line unchanged, quality without its first 6 characters. Read 2 records are unchanged.
- stderr: `Find matched reads: <n>` and `Total reads: <m>`.

**Behaviour notes** — Pairs whose barcode is not in the index file are discarded; there is no "undetermined" output. The barcode length 6 is hard-coded.

**Pitfalls**
- One output file handle is kept open for read 1 and one for read 2 of every barcode at the same time, so the number of barcodes is limited by the open-file limit of the shell.
- The output directory in `<outputPrefix>` must already exist.
- The read loop stops when read 2 ends; file lengths are not compared.
- File names and the prefix are passed to the shell unquoted (`zcat`, `bzcat`, `gzip -f`), so paths with spaces or shell characters break.
- Since the 2026-11-23 fix read 1 keeps its last base; read 1 is 1 bp longer than in output of earlier versions.

**Minimal example** (from `test/unit/perl/run.sh` and the pipeline)
```
perl src/perl/demultiplex_ecker_scWGBS.pl outdir/plate1 barcodes.txt plate1_R1_001.fastq.gz plate1_R2_001.fastq.gz 2> plate1.demultiplex.log
```

### `src/R/call_topdom_on_hic.R`

**Purpose** — Calls topologically associating domains with TopDom on every chromosome of a `.hic` file and writes all domains to one tab-separated file.

**Command** — `Rscript call_topdom_on_hic.R <in.hic> <out.tsv> [chromosome_list]`. It loads `tidyverse`, `TopDom`, `strawr` and `hicrep` with `library()` and calls `logging::loginfo`, so the `logging` package must be installed too. Not verified: this script was not executed in this review and has no unit test.

| Argument | Default | Meaning |
|---|---|---|
| 1 | required | Input `.hic` file. |
| 2 | required | Output file. |
| 3 | `1,2,...,22,X` | Optional comma-separated chromosome names, spelled as in the `.hic` file. An empty string gives the default. |

**Inputs** — A `.hic` file that `hicrep::hic2mat` can read at 50 kb resolution with `KR` normalisation for every listed chromosome (intra-chromosomal matrices only).

**Output** — `<out.tsv>` written with `write_tsv` (tab-separated, with a header line). One row per TopDom domain, all chromosomes together: the columns of TopDom's `bed` table without `name` (the code uses `chrom`, `chromStart`, `chromEnd`), followed by `chrom2`, `chromStart2`, `chromEnd2`, which repeat the first three. Rows that TopDom labels other than `domain` (gaps, boundaries) are removed. Not verified: that the `bed` table has no further columns, that is, that the output has exactly the six columns `chrom, chromStart, chromEnd, chrom2, chromStart2, chromEnd2` (TopDom was not available to check). Log lines are written by `logging`.

**Behaviour notes** — Hard-coded settings: resolution 50,000 bp, normalisation `KR`, TopDom `window.size = 5`, `set.seed(2022)`. For each chromosome the matrix is written to a temporary text file (bin chrom, start, end, then the matrix row; missing values as 0), read with `TopDom::readHiC`, and removed afterwards. The chromosome names in the output are those given in the list.

**Pitfalls**
- The default list has no `chr` prefix; for a `.hic` file with other names, or a non-human genome, pass the third argument.
- Resolution, normalisation and window size cannot be changed from the command line.
- The log line always prints `chr` in front of the name, also when the name already has it.

**Minimal example** (README Analysis 3; the second line shows the optional list from the comment in the script)
```
Rscript src/R/call_topdom_on_hic.R sample.mapQ30.hic sample.mapQ30.topdom_default.tsv
Rscript src/R/call_topdom_on_hic.R in.hic out.tsv 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X
```

---

## Notebooks (not documented in detail)

- `src/python/fit_dist.ipynb` — An exploratory notebook that reads chr20 GCH values from one tabix 6plus2 bed at a hard-coded `/Users/yaping/...` path, z-scores the accessibility (100 minus methylation), fits a five-Gaussian mixture to its histogram with `lmfit`, and writes each site's most likely component to a bedGraph.
- `src/python/predict_tf_binding.ipynb` — An exploratory notebook that trains cross-validated classifiers (scikit-learn and xgboost, with ROC curves) to predict CEBPB binding at motif sites from GCH, HCG, DNase and z-score/binomial feature files at hard-coded `/Users/yaping/...` paths.

---

## Known limitations

Found while reading and running the scripts for this guide (November 2026). They are described here, not fixed.

Hi-C conversion and QC

1. `sam2juicer_new.py` applies no filter and stops with a traceback on an unmapped record; filter the BAM first. Only two records per read name are paired; a third is dropped without a message. Positions are 0-based leftmost.
2. `sam2juicer.py` runs only under Python 2 and never finds contig `MT` in the site file. Use `sam2juicer_new.py`.
3. `mh_reads_summary.v2.py` and `mh_hic_filter.py` count a third record with the same read name as another fragment, and test the duplicate and paired flags of the second record only. BAMs written by `Bhmem` have two records per read name, so this matters only for other aligners.
4. `mh_reads_summary.v2.py` has no unit test; `mh_reads_summary.py` has one. Both expect a spike-in contig named `lambda`.
5. `mh_hic_filter.py` puts every contig outside its exclusion list (`chrUn*`, `*random`, `chrM`, `chrX`, `chrY`) into the per-chromosome test. A short contig such as `lambda` or `chrEBV` with no long cis pair sets `PerChromCisLongFilter` to 1.
6. `--in_cram` is the option name in three scripts, but only BAM input has been exercised.
7. `interact_density_around_regions.py` has default regions at specific chr17 and chr1 loci; `eigen_adjust_direction.py` has default paths that exist only on the original author's machine and a `coverage` parameter that is never set.

GCH/HCG methylation tracks

8. `betabinom_test_mgch.py` imports `fastbetabino` even with `--binomial_only True`, and the fit sample is unseeded for chromosomes with more than 1,000,000 sites. Its comments mention "log2 fold change" and "greater"/"smaller" tests; the code writes fraction differences and offers only `less` and `two-sided`.
9. `smooth_mgch.py`: the comment on `--min_gch` does not match the defaults (100, 1000, 10000). Spline degree 1 and a 2,500 bp gap limit are fixed.
10. `.gz` output of `smooth_mgch.py` and `betabinom_test_mgch.py` is plain gzip, not bgzip, so it cannot be indexed with tabix as written.

Single-cell concordance

11. `cal_concord_sc_at_regions.py` drops a pair without a message when exactly one region ends at or past the chromosome end; `--min_data` means "more than"; all bigWigs are reopened for every BEDPE line. Default permutations: 0.
12. `call_methy_concord_sc.py` changes directory to the input directory, so a relative output path lands there. It skips pairs whose count equals the cell number, and its name column omits the second chromosome. The comment says "at least 3 cells"; the code requires more than 10.
13. `test_indep_poibin.py` ignores the chromosome column when comparing distances and uses a fixed 3 Mb maximum distance. Despite its name it is a tool, not a unit test.

Perl and R

14. `demultiplex_ecker_scWGBS.pl` has no options although the usage shows `[Options]`. The barcode length is fixed at 6. The loop ends on R2 only, compression is detected from the R1 file name, and two file handles per barcode stay open.
15. `call_topdom_on_hic.R` fixes 50 kb resolution, KR normalisation and window 5; the default chromosome list is human 1-22 and X.
