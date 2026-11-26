# Aligners: `Bhmem`, `BhmemMultRef`, `BsmemMultRef`, installation

Part of the [usage guide](USAGE.md). This page is the reference for the three aligners, the installer scripts and the `bisulfitehicMap` wrapper.

- `Bhmem` is the aligner used in practice. It has a README example, a regression test (`test/regression/`) and is what the sc_NOMeHiC pipeline runs.
- `BhmemMultRef` and `BsmemMultRef` (two allele-specific genomes) have no example and no test. Their sections describe what the code does; nothing in them has been run.
- `<repo>` means the checkout directory.
- Sentences that start with "Not verified:" rest on library behaviour that was not checked.

Contents: [`Bhmem`](#bhmem) · [`BhmemMultRef`](#bhmemmultref) · [`BsmemMultRef`](#bsmemmultref) · [Installing and running](#installing-and-running) · [Known limitations](#known-limitations)

---

## `Bhmem`

**Purpose**

`Bhmem` aligns one paired-end bisulfite Hi-C library (Methyl-HiC, NOMe-HiC, single-cell variants) to one reference genome. Each read end is aligned separately with BWA-MEM in single-end mode against the C-to-T and the G-to-A converted genome, and one alignment per end is then chosen as a pair. Choose it whenever there is a single reference; the other two classes are for two allele-specific genomes.

**Command**

```
java -Xmx<N>G \
  -Djava.library.path=<repo>/jbwa/jbwa-1.0.0/src/main/native \
  -cp "<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar:<repo>/jbwa/jbwa-1.0.0/jbwa.jar" \
  main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem \
  <ref.fa> <out.bam> <R1.fq[.gz]> <R2.fq[.gz]> [options]
```

- `-Djava.library.path` must point at the directory that holds `libbwajni.so`; the class calls `System.loadLibrary("bwajni")`.
- `jbwa.jar` must be on the class path in addition to the fat jar. The fat jar in this checkout contains no `com/github/lindenb` classes.
- Positional argument 1, `<ref.fa>`: the unconverted reference FASTA. It is used only to locate the sequence dictionary and the `Bisulfite_Genome` directory; the FASTA itself is not opened.
- Positional argument 2, `<out.bam>`: the complete path of the output BAM. The string is passed unchanged to the BAM writer, so no suffix is added. The usage string calls it `outputPrefix`, which is wrong; the README (`testOut.bam`) and the wrapper (`output.bam`) are right. Exactly one file is written.
- Positional arguments 3 and 4: R1 and R2 FASTQ.
- The usage string puts options before the positional arguments. The README and the regression tests put them after; that order is the one exercised by the tests.
- If fewer than 4 command-line tokens are given (options count as tokens), the program prints the usage line and the option list to stderr and returns normally (exit status 0).

**Inputs**

Reference files, all relative to the directory of `<ref.fa>` (resolved to an absolute path, so a bare file name works):

| File | Used for |
|---|---|
| `Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa` plus its BWA index (`.amb .ann .bwt .pac .sa`) | C-to-T converted genome |
| `Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa` plus its BWA index | G-to-A converted genome |
| a sequence dictionary, first one found of: `<ref without .fa/.fasta/.fna>.dict`, `<ref.fa>.dict`, `<ref.fa>.fai` | `@SQ` lines of the BAM header |

- The file names under `Bisulfite_Genome` are fixed in the code. They are the names `bismark_genome_preparation` produces (README); each converted FASTA must then be indexed with `bwa index`.
- Contig names in the converted genomes may carry the suffix `_CT_converted` or `_GA_converted`; the suffix is removed from every alignment. After removal the name must exist in the sequence dictionary.
- If no dictionary file is found the program stops with an `IOException` that names the missing files.

FASTQ:

- Files are opened with zlib `gzopen`, so plain and gzip-compressed FASTQ both work. The file name `-` means standard input.
- R1 and R2 are read in step, one record from each. Reading stops when either file ends; extra records in the longer file are ignored without a message.
- The read name is the FASTQ name up to the first white space. R1 and R2 alignments are joined by read name, so the two names must be identical. BWA is called through the JNI layer, which does not strip `/1` and `/2`.
- `-readNameSuffix <regex>`: a Java regular expression; every match is deleted from the read name of each alignment (`String.replaceAll`). Example values from the source comment: `_\d$`, `_R\d$`. The pattern is not anchored unless you anchor it.

`-enzymeList` file:

- Tab-separated text, optionally gzip-compressed (detected by the `.gz` suffix only).
- Column 1 contig, column 2 start, column 3 end, read as BED (0-based, half-open). Further columns are ignored. Lines that start with `#` are skipped.
- Any other non-data line (blank line, `track` line) causes a parse exception.
- Contig names must be the reference contig names (without the `_CT_converted` suffix); lookup is case-sensitive.

**Options**

| Option | Default | What the code does |
|---|---|---|
| `-clip5` | 0 | Sets BWA-MEM's 5'-end clipping penalty (`pen_clip5`, the first value of `bwa mem -L`). It does not trim bases, contrary to the usage text "bp at 5' end to clip". The value is always passed, so the default 0 replaces BWA's own default of 5. |
| `-clip3` | 0 | Same for the 3'-end clipping penalty (`pen_clip3`). |
| `-buffer` | 100000 | Batch size in read pairs. A batch is aligned when it holds `buffer + 1` pairs or at end of input. Also the period of the progress message, there counted in reads (two per pair). `-buffer 0` causes a division by zero. |
| `-rgId` | null | With `-rgSm`: adds `@RG ID:<rgId> LB:<rgId> SM:<rgSm> PL:illumina PU:barcode`. Independently of `-rgSm`, every record gets `RG:Z:<rgId>` when `-rgId` is given. |
| `-rgSm` | null | Sample name for the `@RG` line. The `@RG` line is written only when both `-rgId` and `-rgSm` are non-empty. |
| `-noSoftClipOnSuppl` | false | When absent, BWA's `MEM_F_SOFTCLIP` flag is set (`bwa mem -Y`): supplementary alignments are soft-clipped. When present the flag is not set and BWA hard-clips supplementary alignments; SEQ and QUAL of such a record then contain only the aligned part. |
| `-nonDirectional` | false | Adds two more kinds of candidate pair: R1 from the CT genome with R2 from the GA genome, and R1 from the GA genome with R2 from the CT genome. It does not add any in-silico conversion or any alignment run. |
| `-forceSameChr` | false | Candidate pairs are formed only from an R1 alignment and an R2 alignment on the same contig. Reads with no such combination are dropped. |
| `-outputMateDiffChr` | false | When absent, a read whose chosen pair lies on two different contigs is dropped (both records). When present it is written. |
| `-enzymeList` | null | Restriction-site intervals used as a tie-break between candidate pairs whose MAPQ sums are equal (step 6 below). |
| `-readNameSuffix` | null | Regular expression deleted from read names (see Inputs). |
| `-pbat` | false | R1 is converted G-to-A and R2 C-to-T (the reverse of the default), and the strand flag (0x10) of every record is inverted. SEQ and QUAL are not changed by this option. |
| `-score` | 0 | BWA-MEM's minimum score to output (`opt->T`, `bwa mem -T`). The default 0 replaces BWA's own default of 30. |
| `-t` | 1 | BWA worker threads (`opt->n_threads`) used inside each alignment call. |
| `-h` | false | No effect. The flag is tested before the command line is parsed, so it is always false at that point, and it is not read afterwards. The option list is printed only when fewer than 4 tokens are given. |

**How it aligns**

1. Read pairs are read into a batch of `buffer + 1` pairs.
2. In-silico conversion, one per read end:

   | Mode | R1 | R2 |
   |---|---|---|
   | default | C to T | G to A |
   | `-pbat` | G to A | C to T |
   | `-nonDirectional` | no change to the conversions | no change |

3. Each converted read is aligned in single-end mode to both converted genomes, in every mode. That is four BWA calls per batch (R1 on CT, R1 on GA, R2 on CT, R2 on GA), run one after another. BWA parameters set by this program: threads (`-t`), minimum output score (`-score`), soft-clipping of supplementary alignments (on unless `-noSoftClipOnSuppl`), clipping penalties (`-clip5`, `-clip3`). All other BWA parameters are the `mem_opt_init` defaults.
4. Every SAM line BWA returns becomes a candidate record, except lines flagged unmapped, secondary (0x100) or QC-fail. Supplementary lines (0x800) are kept as candidates and their supplementary flag is cleared. A read that spans a ligation junction therefore contributes its primary part and its supplementary part as separate candidates for that end.
5. Candidate pairs for a read are all combinations of one R1 candidate and one R2 candidate from the same converted genome (CT with CT, GA with GA). With `-nonDirectional` the cross-genome combinations are added. With `-forceSameChr` only combinations on the same contig are formed. There is no strand requirement on the pair.
6. The first candidate pair becomes the current best. Each further candidate (the challenger) replaces the current best when the first of these tests that applies says so:
   1. The challenger has MAPQ > 0 on both ends and the current best has MAPQ 0 on at least one end: challenger wins.
   2. Otherwise the larger MAPQ sum (R1 + R2) wins.
   3. MAPQ sums equal: if the current best is on two contigs and the challenger is on one contig, the challenger wins. (This test is applied in that direction only.)
   4. `-enzymeList` given: if exactly one of the two pairs has both ends on contigs that occur in the list, that pair wins. Otherwise, if exactly one of the two pairs has both ends within 50 bp of a listed site (alignment start - 50 to alignment end + 50 overlaps a site), that pair wins.
   5. Larger sum of `AS` wins.
   6. Smaller sum of `NM` wins.
   7. Larger total number of aligned `M` bases in the two CIGARs wins.
   8. Still equal: the current best is kept.
7. The chosen pair is the only output for the read: one record for R1 and one for R2. The parts of a chimeric read that were not chosen are not written; they remain visible only as the soft-clipped part of SEQ (or are absent when the chosen record was hard-clipped under `-noSoftClipOnSuppl`). No `SA` or `XA` tag is kept.
8. If the two records are on different contigs and `-outputMateDiffChr` is absent, the pair is discarded.
9. Mate fields, proper-pair flag and TLEN are set (see Output), and the pair is written.

**Output**

One BAM file at the path given as argument 2.

- Sort order: header says `@HD VN:1.5 SO:unsorted`. Records are written batch by batch; inside a batch the order of reads is the iteration order of a hash map, not the input order. The two records of a pair are adjacent, R1 first.
- Header: `@SQ` lines from the sequence dictionary; one `@RG` line when both `-rgId` and `-rgSm` are given; no `@PG` line.
- Records written: only reads for which both ends have at least one usable alignment and a candidate pair exists. Dropped without a record: reads with either end unmapped; reads with no same-genome combination (unless `-nonDirectional`); reads with no same-contig combination under `-forceSameChr`; pairs on different contigs unless `-outputMateDiffChr`. There are no unmapped, secondary or supplementary records in the output, and no MAPQ filter: pairs with MAPQ 0 are written.
- Flags: 0x1 always; 0x40 on R1, 0x80 on R2; 0x10 as described below; 0x20 copies the mate's 0x10; 0x8 never; 0x800 and 0x100 never; 0x400 never (no duplicate marking). 0x2 (proper pair) is set when both records are on the same contig and their strand flags differ (orientations FR and RF are both accepted); there is no distance limit. Same-strand pairs and inter-contig pairs have 0x2 clear.
- RNEXT, PNEXT and the tags `MC` and `MQ` are filled from the mate.
- Strand flag, default mode: 0x10 is what BWA reported for the converted read against the converted genome. SEQ is the original (unconverted) read when 0x10 is clear and its reverse complement when 0x10 is set.
- Strand flag, `-pbat`: 0x10 is the inverse of what BWA reported. SEQ and QUAL are still stored on the reference forward strand, so a record with 0x10 set holds the read as it is in the FASTQ and a record with 0x10 clear holds its reverse complement.
- SEQ: the unconverted read from the FASTQ, on the reference forward strand, with hard-clipped bases removed. QUAL: taken unchanged from BWA's output (BWA already reversed it for reverse-strand alignments and trimmed it for hard clips).
- MAPQ: the value BWA-MEM computed for the single-end alignment of the converted read to the one converted genome the record comes from. It is not recomputed across the two converted genomes. For a record that was a supplementary alignment, BWA caps the value at the primary alignment's MAPQ.
- TLEN: first the htsjdk value: 0 if the contigs differ, otherwise the signed distance between the two 5' ends, inclusive (5' end = POS for a record with 0x10 clear, alignment end for a record with 0x10 set; soft clips not counted). Then TLEN is set to 0 on both records if any of these holds: different contigs; equal strand flags; the record with 0x10 set starts at a lower POS than the record with 0x10 clear. TLEN is therefore non-zero only for same-contig pairs with opposite strand flags in which the 0x10-clear record starts at or before the 0x10-set record. The test uses the flags as written, that is after the `-pbat` inversion (see Known limitations).
- Tags:

  | Tag | Meaning |
  |---|---|
  | `XG:Z` | `CT` or `GA`: the converted genome the record was aligned to. |
  | `XR:Z` | `CT` or `GA`: the in-silico conversion applied to the read. Default mode: R1 `CT`, R2 `GA`. `-pbat`: R1 `GA`, R2 `CT`. |
  | `AS:i` | BWA alignment score, copied. |
  | `NM:i` | BWA edit distance of the converted read against the converted genome, copied. It does not count bisulfite mismatches of SEQ against the unconverted reference. |
  | `MD:Z` | BWA's MD string, also in converted space, copied. |
  | `XS:i` | BWA suboptimal alignment score, copied when BWA wrote it. |
  | `RG:Z` | Value of `-rgId`, when given. |
  | `MC:Z` | Mate CIGAR, added by htsjdk. |
  | `MQ:i` | Mate MAPQ, added by htsjdk. |

  Every other tag BWA writes is copied too, except `SA` and `XA`, which are removed. Integer tags stay integers; all other types are stored as strings.
- Log (log4j console output, lines of the form `<ms> [main] INFO <class>  - <message>`):
  - `Parsing input -enzymeList bed file ...` when the option is given.
  - `Processing reads <N> ...` and `Unique mapped and paired reads <M> ...` each time the number of reads read is a multiple of `-buffer`.
  - At the end: `<N> reads in total`; `Unique Mapped and paired <M> reads (<pct>%)`; `Bhmem's running time is: ...`.
  - `N` counts reads (2 per pair read). `M` counts records written (2 per pair written). No uniqueness or MAPQ test is behind the word "Unique".
  - BWA's own messages (`[M::bwa_idx_load_from_disk]`, `[M::mem_process_seqs] Processed ...`) go to stderr.
  - Not verified: the log4j lines go to stdout (log4j `BasicConfigurator` default; the log4j source was not read, and the regression runner merges both streams).

**Pitfalls**

- The second argument is a file name, not a prefix. An existing file at that path is overwritten.
- The output is unsorted and has no duplicate marks. The README pipeline before methylation calling is `samtools sort -n | samtools fixmate -m | samtools sort | samtools markdup | samtools calmd` followed by `samtools index`. `calmd` matters because `NM` and `MD` are in converted space.
- Read names of R1 and R2 must match exactly. With `/1` and `/2` suffixes and no `-readNameSuffix`, no pair is formed and the BAM is empty.
- `-clip5` and `-clip3` are penalties, not trimming lengths. Leaving them at 0 is not the BWA default.
- `-score 0` (the default) makes BWA report alignments that `bwa mem` would suppress with its default `-T 30`.
- Memory: two BWA indexes are loaded at the same time by native code (`bwa_idx_load`), outside the Java heap set with `-Xmx`. The Java heap holds one batch: `buffer + 1` read pairs, four arrays of SAM text and one record object per candidate line. The README uses `-Xmx20G` with `-buffer 100000` and `-Xmx15G` with `-buffer 1000`.
- `-t` only parallelises the inside of each BWA call. The four calls per batch are sequential and the pair selection is single-threaded.
- Without `-outputMateDiffChr` all inter-chromosomal contacts are discarded.
- Pairs with MAPQ 0 are in the output; filter downstream.
- TLEN is 0 for most Hi-C pairs. Use POS and PNEXT for distances.
- `-rgId` without `-rgSm` writes `RG` tags that refer to an `@RG` line that is not in the header.
- A usage error (too few arguments, unknown option) prints the usage and exits with status 0.
- Directional mode converts R1 C-to-T and R2 G-to-A only. `-nonDirectional` does not add the opposite conversions; `-pbat` swaps them.

---

## `BhmemMultRef`

**Purpose**

`BhmemMultRef` aligns one paired-end bisulfite library to two allele-specific genomes (for example two parental strains), aligning each read end separately as `Bhmem` does, and assigns each read pair to the genome where it aligns better. It keeps only pairs whose two ends are on the same contig with opposite strands. Choose it over `Bhmem` when two personal genomes with chain files to a common reference are available; `BsmemMultRef` is the variant that uses BWA's paired-end mode.

**Command**

```
java -Xmx<N>G \
  -Djava.library.path=<repo>/jbwa/jbwa-1.0.0/src/main/native \
  -cp "<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar:<repo>/jbwa/jbwa-1.0.0/jbwa.jar" \
  main.java.edu.mit.compbio.bisulfitehic.mapping.BhmemMultRef \
  <ref.fa> <outputPrefix> <R1.fq[.gz]> <R2.fq[.gz]> \
  -refs <dir/genome1.fa> -refs <dir/genome2.fa> \
  -chains <genome1_to_ref.chain> -chains <genome2_to_ref.chain> [options]
```

This command shape is built from the source; the README has no example and no test runs this class.

- Positional argument 1, `<ref.fa>`: the common reference. Only its sequence dictionary is used (header of the merged and ambiguous BAMs). No BWA index of `<ref.fa>` is loaded.
- Positional argument 2, `<outputPrefix>`: a real prefix. Four files are written: `<outputPrefix>.merged.bam`, `<outputPrefix>.ambiguous.bam`, `<outputPrefix>.alignToRef_1.bam`, `<outputPrefix>.alignToRef_2.bam`. The number in `alignToRef_N` is the position of the genome in the `-refs` list.
- Positional arguments 3 and 4: R1 and R2 FASTQ.
- `-refs` and `-chains` are required, exactly two of each, repeated as separate options (form taken from the usage string). The n-th `-chains` file belongs to the n-th `-refs` genome.
- Fewer than 3 tokens prints the usage to stderr and returns (exit status 0). A missing or mismatched `-refs`/`-chains` list is not handled that way: it raises an uncaught `CmdLineException`. A `-refs` count other than two raises `IllegalArgumentException("Number of reference genome is not equal to two!!")`.

**Inputs**

- For each `-refs` genome: a sequence dictionary found by the same rule as in `Bhmem` (`genome.dict`, `genome.fa.dict` or `genome.fa.fai`), and `<parent of the -refs path>/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa` and `.../GA_conversion/genome_mfa.GA_conversion.fa` with BWA indexes. The parent directory is taken from the path as typed, so each `-refs` value must contain a directory component; a bare file name produces the path `null/Bisulfite_Genome/...`. The two genomes must be in different directories, because the converted-genome file names are fixed.
- For `<ref.fa>`: a sequence dictionary only.
- `-chains`: chain files read with htsjdk `LiftOver`, used to lift positions from each allele genome to `<ref.fa>` coordinates, with a minimum match of 0.9.
- FASTQ: same reader as `Bhmem` (plain or gzip, `-` for stdin, reading stops at the shorter file). Read names of R1 and R2 must be identical. `-trimReadName` removes a trailing `/<digits>` from both names before alignment. There is no `-readNameSuffix` and no `-enzymeList`.

**Options**

| Option | Default | What the code does |
|---|---|---|
| `-refs` | null | Allele genome FASTA; give the option twice. Exactly two are required. |
| `-chains` | null | Chain file for the genome at the same position in `-refs`; give the option twice. Required even with `-addAmbiguousReads`, where the chains are loaded but not used. |
| `-clip5` | 0 | BWA 5' clipping penalty, as in `Bhmem`; not a trimming length. |
| `-clip3` | 0 | BWA 3' clipping penalty, as in `Bhmem`. |
| `-buffer` | 100000 | Batch size in read pairs (`buffer + 1` per batch) and period of the progress message. |
| `-rgId` | null | As in `Bhmem`; applied to all four BAM headers and to every record. |
| `-rgSm` | null | As in `Bhmem`. |
| `-noHicMode` | false | No effect. The field is never read. The usage text itself says "not implemented yet". |
| `-addAmbiguousReads` | false | Switches to a different code path: for each genome, every joined pair is written to that genome's BAM. No best-genome choice is made, nothing is written to the merged and ambiguous BAMs (they contain only a header), and `-ratio`, `-seed` and the chain files are not used. |
| `-noChimericMode` | false | No effect. The field is never read. The usage text says "not implemented yet". |
| `-trimReadName` | false | Deletes the pattern `/\d+$` from R1 and R2 read names before alignment. |
| `-ratio` | 0 | Threshold on the difference `(MAPQ R1 + MAPQ R2 of the best genome) - (MAPQ R1 + MAPQ R2 of the second genome)`. The pair is written to the best genome's BAM only if the difference is strictly greater than `-ratio`. It is not a probability ratio. |
| `-score` | 0 | BWA minimum output score, as in `Bhmem`. |
| `-seed` | -1 | Seed of the Mersenne Twister used for the coin flip in the genome comparison; -1 means unseeded. |
| `-t` | 1 | BWA worker threads per alignment call. |
| `-h` | false | No effect, for the same reason as in `Bhmem`. |

**How it aligns**

1. Read pairs are batched as in `Bhmem`. R1 is converted C-to-T and R2 G-to-A. There is no `-pbat` and no `-nonDirectional` mode.
2. For each of the two genomes, each converted read is aligned in single-end mode to that genome's CT and GA converted genomes: eight BWA calls per batch. Soft-clipping of supplementary alignments is always on.
3. Within one genome, candidate lines are filtered as in `Bhmem` (unmapped, secondary and QC-fail removed; supplementary kept with the flag cleared). For each combination of read, contig, read end and converted genome, one alignment is kept: higher MAPQ; then more `M` bases; then lower `NM`; then higher `AS`; if still equal, the one seen first.
4. A candidate pair is an R1 and an R2 alignment on the same contig, from the same converted genome, with opposite strands. Inter-contig pairs and same-strand pairs are never formed. The best pair per read is chosen by: higher minimum MAPQ of the two ends; then higher MAPQ sum; then more `M` bases; then lower `NM` sum; then higher `AS` sum; if still equal, the one seen first.
5. Default path (no `-addAmbiguousReads`): the best pair from genome 1 is the initial best. The pair from genome 2 is compared with it on, in order, minimum MAPQ, MAPQ sum, `M` bases and `NM` sum.
   - Genome 2 better on the first criterion that differs: genome 2 becomes best, genome 1 second.
   - The four criteria equal and genome 2 has the higher `AS` sum: a coin flip. On one outcome genome 2 becomes best and genome 1 second; on the other outcome the genome 2 pair is forgotten and no second-best is recorded.
   - Otherwise (genome 1 better, or all five criteria equal): genome 1 stays best and genome 2 becomes second.
6. The best pair is written to `alignToRef_<best>.bam` if no second-best pair is recorded, or if the MAPQ-sum difference between best and second is greater than `-ratio`.
7. Independently of step 6, the alignment start of each end of the best pair is lifted to `<ref.fa>` coordinates with the best genome's chain. If either start cannot be lifted the pair is not written to the merged or ambiguous BAM and is not counted as mapped. Otherwise POS is replaced by the lifted position, mate fields, proper-pair flag and TLEN are recomputed, and the pair is written to both `merged.bam` and `ambiguous.bam`. The contig name and the CIGAR are not changed by the lift.
8. `-addAmbiguousReads` path: steps 5 to 7 are replaced by writing every pair from step 4 to the BAM of the genome it was aligned to.

**Output**

- Four BAM files (see Command), all written with the same writer settings as `Bhmem`: unsorted, pairs adjacent with R1 first, batch order with hash-map order inside a batch, no `@PG` line, `@RG` as in `Bhmem`.
- `alignToRef_N.bam`: header from the dictionary of the N-th `-refs` genome; coordinates in that genome. Contains the pairs assigned to genome N by step 6 (or all pairs of genome N with `-addAmbiguousReads`).
- `merged.bam`: header from the dictionary of `<ref.fa>`; the best pair of every read that could be lifted, whether or not it was assigned to a genome in step 6.
- `ambiguous.bam`: receives exactly the same records as `merged.bam`. The code has no separate test for "ambiguous".
- Records: both ends mapped, same contig, opposite strands. Reads failing that are dropped. No unmapped, secondary or supplementary records.
- Flags: 0x1, 0x40/0x80, 0x10 as BWA reported, 0x20 from the mate, 0x2 by the same rule as `Bhmem` (same contig and orientation FR or RF). SEQ and QUAL as in `Bhmem` default mode.
- MAPQ: BWA single-end MAPQ in the genome the record comes from.
- TLEN: the htsjdk 5'-to-5' value. There is no chimeric reset in this class.
- Tags: BWA's tags except `SA` and `XA` (`NM`, `MD`, `AS`, `XS`, in converted space), `RG` when `-rgId` is given, and `MC`, `MQ` from htsjdk. `XG` and `XR` are not written. No tag records the assigned genome; the assignment is expressed only by which file a pair is in.
- Log: progress lines as in `Bhmem`, then at the end `<N> reads in total` and six counts with percentages of total reads: `Unique Mapped and paired`; `... paired in the same chr`; `... same chr and MapQ >=1`; `... MapQ >=30`; `... MapQ >=1 and insertionSize > 20000`; `... MapQ >=30 and insertionSize > 20000`; then `BhmemMultFa's running time is: ...`. In the default path the counts are taken from the pairs written to `merged.bam`; with `-addAmbiguousReads` they are taken from genome 1 only. "MapQ >=" tests the smaller MAPQ of the two ends; "insertionSize" is `|TLEN|` of R1.
- Not verified: whether records in `merged.bam` get their reference index re-resolved against the `<ref.fa>` dictionary depends on htsjdk `SAMRecord.setHeader`, which was not read. Using allele genomes whose contig names and order equal those of `<ref.fa>` avoids the question.

**Pitfalls**

- No README example and no test covers this class.
- Each `-refs` path needs a directory component, and each genome needs its own directory with its own `Bisulfite_Genome`.
- Four BWA indexes are loaded at once by native code, outside the Java heap.
- `ambiguous.bam` is a copy of `merged.bam`, not a list of ambiguous reads. To get the ambiguous reads, take the reads in `merged.bam` that are in neither `alignToRef` file.
- With the default `-ratio 0`, a pair that maps to both genomes with equal MAPQ sums is written to neither `alignToRef` file.
- Inter-chromosomal pairs are never reported.
- Output is unsorted and unmarked for duplicates, as with `Bhmem`.
- `-clip5`, `-clip3`, `-score`, `-buffer`, `-t`, memory and exit-status remarks from `Bhmem` apply.

---

## `BsmemMultRef`

**Purpose**

`BsmemMultRef` aligns one paired-end bisulfite library to two allele-specific genomes using BWA-MEM's paired-end mode (R1 and R2 aligned together) and assigns each pair to the better genome with the same logic as `BhmemMultRef`. The `-noHicMode` usage text of `BhmemMultRef` names this class as the WGBS mode. Choose it for two-genome alignment of ordinary paired-end bisulfite libraries; choose `BhmemMultRef` when ends must be aligned independently.

**Command**

```
java -Xmx<N>G \
  -Djava.library.path=<repo>/jbwa/jbwa-1.0.0/src/main/native \
  -cp "<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar:<repo>/jbwa/jbwa-1.0.0/jbwa.jar" \
  main.java.edu.mit.compbio.bisulfitehic.mapping.BsmemMultRef \
  <ref.fa> <outputPrefix> <R1.fq[.gz]> <R2.fq[.gz]> \
  -refs <dir/genome1.fa> -refs <dir/genome2.fa> \
  -chains <genome1_to_ref.chain> -chains <genome2_to_ref.chain> [options]
```

Positional arguments, required options, output file names (`<outputPrefix>.merged.bam`, `.ambiguous.bam`, `.alignToRef_1.bam`, `.alignToRef_2.bam`) and error behaviour are the same as for `BhmemMultRef`. The usage string printed by this class says `BhmemMultRef`. The command shape is built from the source; no README example and no test runs this class.

**Inputs**

Same as `BhmemMultRef`: dictionaries and `Bisulfite_Genome` BWA indexes for both `-refs` genomes (with the same directory-component requirement), a dictionary for `<ref.fa>`, two chain files, FASTQ plain or gzip. R1 and R2 are paired by position in the files, not by name; the R1 name is the key of the pair. `-trimReadName` removes a trailing `/<digits>` from both names.

**Options**

| Option | Default | What the code does |
|---|---|---|
| `-refs` | null | As in `BhmemMultRef`: two required. |
| `-chains` | null | As in `BhmemMultRef`: two required, also with `-addAmbiguousReads`. |
| `-clip5` | 0 | BWA 5' clipping penalty; not a trimming length. |
| `-clip3` | 0 | BWA 3' clipping penalty. |
| `-buffer` | 100000 | Batch size in read pairs (`buffer + 1`) and period of the progress message. Each batch is one BWA paired-end call per converted genome, and BWA estimates the insert-size distribution from that batch alone. |
| `-rgId` | null | As in `Bhmem`. |
| `-rgSm` | null | As in `Bhmem`. |
| `-chimericMode` | false | No effect. The field is never read. The usage text says "not implemented yet". |
| `-addAmbiguousReads` | false | As in `BhmemMultRef`: every pair to every genome's BAM, nothing to merged/ambiguous, no use of `-ratio`, `-seed` or the chains. |
| `-trimReadName` | false | Deletes `/\d+$` from read names before alignment. |
| `-ratio` | 0 | MAPQ-sum difference threshold, as in `BhmemMultRef`. |
| `-score` | 0 | BWA minimum output score. |
| `-seed` | -1 | Seed for the coin flip, as in `BhmemMultRef`. |
| `-t` | 1 | BWA worker threads per alignment call. |
| `-h` | false | No effect. |

**How it aligns**

1. Batching and conversion as in `BhmemMultRef`: R1 C-to-T, R2 G-to-A.
2. For each genome, the batch is aligned in BWA paired-end mode (`MEM_F_PE` set) to the CT genome and to the GA genome: four BWA calls per batch. Soft-clipping of supplementary alignments is always on. No insert-size distribution is passed in, so BWA infers it per call.
3. For each read end only the first SAM line BWA returns is used; further lines (for example supplementary alignments) are ignored. The pair is skipped for that converted genome if either first line is unmapped, secondary or QC-fail.
4. Mate fields and the proper-pair flag are set; then the CT-genome pair and the GA-genome pair of the same read are compared: higher minimum MAPQ; then higher MAPQ sum; then more `M` bases; then lower `NM` sum; then higher `AS` sum; if still equal, the CT-genome pair is kept. There is no same-contig or opposite-strand requirement in this class.
5. Genome comparison, assignment to `alignToRef_N.bam`, lift-over, and writing of `merged.bam` and `ambiguous.bam` are identical to steps 5 to 8 of `BhmemMultRef`.

**Output**

Same four files, header rules, ordering, tags (no `XG`/`XR`), MAPQ source and log lines as `BhmemMultRef`, with these differences:

- A pair may lie on two contigs or have both ends on the same strand; such pairs have the proper-pair flag clear and, for different contigs, TLEN 0.
- MAPQ and the BWA tags come from BWA's paired-end alignment.
- The log line "paired in the same chr" is a real subset here.

**Pitfalls**

- All pitfalls of `BhmemMultRef` except the one about inter-chromosomal pairs.
- Supplementary alignments are discarded, so this class does not recover split (chimeric) reads.
- A small `-buffer` gives BWA few pairs to estimate the insert-size distribution from.

---

## Installing and running

**What `install.sh` does**

1. `cd jbwa/jbwa-1.0.0/ && make && cd ../../ && mvn compile assembly:single`, as one `&&` chain.
   - `make` (default target `all`) in `jbwa/jbwa-1.0.0`: aborts if `JAVA_HOME` is not defined; downloads the BWA source with `wget` from GitHub and builds it if `bwa-9f26bfcc.../libbwa.a` is missing; compiles the jbwa Java classes with `javac`; generates `src/main/native/bwajni.h` with `javah`; compiles `bwajni.c` and the BWA objects with `-fPIC` using `$CC` (default `gcc`) and the JNI headers under `$JAVA_HOME/include`; links `src/main/native/libbwajni.so` against zlib and pthread; indexes `test/ref.fa` and runs two small test programs; creates `jbwa.jar` and `jbwa-native.tar`.
   - `mvn compile assembly:single` in the repository root builds `target/bisulfitehic-0.38-jar-with-dependencies.jar` (artifact `bisulfitehic`, version `0.38`, manifest main class `...mapping.Bhmem`). Maven dependencies: args4j 2.0.31, htsjdk 2.13.1, commons-math3 3.3, log4j 1.2.17, commons-lang3 3.4, commons-collections4 4.1, biojava-genome 4.1.0, picard 1.141, and local jars with system scope (`jbwa/jbwa-1.0.0/jbwa.jar`, `lib/java-genomics-io.jar`, `lib/ccInference-0.19-slim.jar`, `lib/gatk-package-distribution-3.3.jar`, `lib/JavaMI.jar`). `jbwa.jar` must exist before Maven runs, which is why `make` comes first.
2. If the environment variable `BISHIC` is empty: `export BISHIC=${PWD}`, append that export line to `~/.bash_profile`, and source `~/.bash_profile`. If `BISHIC` is already set it is left alone.
3. `chmod 755 bisulfitehicMap`.
4. `cd ./target && ln -s $(ls -t | head -n1) bisulfitehic-default.jar`: a symlink named `bisulfitehic-default.jar` to whatever entry of `target/` was modified most recently.

The script has no `set -e`; later lines run even if the first chain fails. The README runs it with `source ./install.sh`, which leaves the calling shell in `target/`.

`install.sudo.sh` is the same with `sudo mvn` and `sudo ln`, in a different order: it changes into `target/` to make the symlink before the `BISHIC` block, so when `BISHIC` is empty it is set to `<repo>/target`, not to `<repo>`.

**Environment the scripts ask for**

| Source | Requirement |
|---|---|
| README | Oracle JDK 8, GNU Make 3.81, gcc 4.8.2, wget, maven, python 3.6, `bismark_genome_preparation` |
| jbwa `Makefile` | `JAVA_HOME` defined and containing `include/*.h`; `javac`, `javah`, `java`, `jar` on `PATH`; `gcc`; `wget` and `unzip` when BWA has to be downloaded; zlib |
| `scripts/quest_install.sh` | modules `java/jdk1.8.0_191`, `apache-maven-mvn3/3.9.7`, `gcc/4.8.3` (and `samtools/1.10.1`, runtime only); `git`; the `module` command; a JDK that still has `javah` (the script checks for it and states that JDK 8 must be first on `PATH`) |
| `test/regression/run_tests.sh` | `java` (loads module `java/jdk-17.0.2+8` if none on `PATH`), `samtools`, `python3`, the bundled `bwa` and `libbwajni.so` |

The jar currently in `target/` records `Build-Jdk-Spec: 11` in its manifest.

**What `scripts/quest_install.sh` adds**

It is a full reinstall for the Quest cluster, configured by environment variables (`TARGET`, `REPO_URL`, `REPO_REF`, `MOD_JAVA8`, `MOD_MAVEN`, `MOD_GCC`, `MOD_SAMTOOLS`, `TAG_JAR`, `JAR_VERSION_TAG`, `RUN_TESTS`). The default `TARGET` is set at the top of the script. Steps:

1. Pre-flight: parent directory exists and is writable; an existing `TARGET` must contain `install.sh`; today's backup directory must not exist; `git` and `module` are available; the three build modules exist.
2. Moves an existing `TARGET` to `${TARGET}_old_<YYYYMMDD>`.
3. `git clone` of `REPO_URL` into `TARGET`, and `git checkout REPO_REF` if it is not `master`. This needs network access (the script says to run it on a login node).
4. `module purge`, then loads Maven, gcc, and JDK 8 last; sets `JAVA_HOME` from the `java` on `PATH`; checks `jni.h`, `javah`, and that `javac` belongs to the same JDK.
5. `bash install.sh`.
6. Checks that `target/bisulfitehic-0.38-jar-with-dependencies.jar` exists; that `Bhmem.class` contains `isLikelyChimericLigation` and the constants `XG` and `XR`; that `lib/ccInference-0.19-slim.jar` exists and is under 1 MB.
7. Runs `test/regression/run_tests.sh` (log in `test/regression/last_run.log`) unless `RUN_TESTS=0`.
8. Copies the jar to `target/bisulfitehic-0.38-NOMeHiC_v1.0-jar-with-dependencies.jar` unless `TAG_JAR=0`.
9. Prints a summary with a rollback command.

Because step 5 uses `bash install.sh`, the `export BISHIC` inside `install.sh` does not reach the calling shell; only the line appended to `~/.bash_profile` persists, and only if `BISHIC` was empty.

**Where things end up**

| Item | Path |
|---|---|
| Fat jar | `<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar` |
| Symlink for the wrapper | `<repo>/target/bisulfitehic-default.jar` |
| JNI library | `<repo>/jbwa/jbwa-1.0.0/src/main/native/libbwajni.so` |
| jbwa classes | `<repo>/jbwa/jbwa-1.0.0/jbwa.jar` |
| Bundled `bwa` binary (for `bwa index`) | `<repo>/jbwa/jbwa-1.0.0/bwa-9f26bfcc7780753129b60717ecab0ebba6f04b7c/bwa` |

**The `bisulfitehicMap` wrapper**

A Perl script that builds and runs one `Bhmem` command:

```
bisulfitehicMap [Options] genome.fa output.bam r1.fq[.gz] r2.fq[.gz]
```

- It reads the install directory from the environment variable `BISHIC`.
- It runs `java -Xmx<mem>G -Djava.library.path=$BISHIC/jbwa/jbwa-1.0.0/src/main/native/ -cp "$BISHIC/target/bisulfitehic-default.jar:$BISHIC/jbwa/jbwa-1.0.0/jbwa.jar" main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem genome.fa output.bam r1 r2 -t <cpu> -score <score> -buffer <buffer> -clip5 <clip5> -clip3 <clip3>`, and adds `-rgId <rgId> -rgSm <rgSm>` only when both are given.
- Options it exposes:

  | Wrapper option | Default | Passed as |
  |---|---|---|
  | `--cpu` | 1 | `-t` |
  | `--mem` | 10 | `-Xmx<mem>G` |
  | `--buffer` | 100000 | `-buffer` |
  | `--score` | 0 | `-score` |
  | `--clip5` | 0 | `-clip5` |
  | `--clip3` | 0 | `-clip3` |
  | `--rgId` | empty | `-rgId` |
  | `--rgSm` | empty | `-rgSm` |

- It has no way to pass `-pbat`, `-nonDirectional`, `-outputMateDiffChr`, `-forceSameChr`, `-enzymeList`, `-readNameSuffix` or `-noSoftClipOnSuppl`. Only the first four positional arguments are used.
- It prints the command to stderr and dies if the command returns non-zero.

Would it work as written:

- It needs `BISHIC` to be the repository root. If `BISHIC` is unset the paths start at `/jbwa/...` and `/target/...`. After `install.sudo.sh` on a machine where `BISHIC` was empty, `BISHIC` is `<repo>/target`, which is wrong for the wrapper.
- It needs `target/bisulfitehic-default.jar` to point at the fat jar. That depends on which entry of `target/` was newest when `install.sh` ran. On the install this guide was written from, the symlink points at the directory `target/classes`, not at the jar. `ln -s` without `-f` does not replace an existing symlink, so re-running `install.sh` does not correct it.
- Not verified: with the symlink pointing at `target/classes`, the wrapper is expected to fail with a missing-class error, because that directory holds the project's own classes but not htsjdk, args4j or log4j. The wrapper was not run.
- Calling `java` directly with the full jar name, as in the README and in `test/regression/run_tests.sh`, does not depend on the symlink or on `BISHIC`.

---

## Known limitations

Found while reading the source for this guide (November 2026). They are described here, not fixed.

`Bhmem`

1. **Strand flags, proper-pair flag and TLEN with `-pbat`.** `-pbat` inverts the strand flag of every record so that methylation callers read the right bisulfite strand. Anything that derives pair orientation from the flags then sees it mirrored: a pair whose reads face each other looks outward-facing (RF), and the reverse.
   - In the raw BAM, the chimeric TLEN reset is applied to the inverted flags, so ordinary inward-facing pairs get TLEN 0 and outward-facing pairs keep it. In the regression BAMs all 100 simulated ordinary PBAT pairs have TLEN 0.
   - `samtools fixmate` recomputes TLEN for every same-contig pair and clears the proper-pair flag of every pair that is not FR by its flags (checked on the regression BAMs). After the standard post-processing, the proper-pair flag therefore no longer means "plausible pair". In one production scNOMe-HiC cell (non-duplicate, MAPQ at least 30), 55% of pairs had no proper-pair flag. Most of those (44% of all pairs) were pairs less than 1 kb apart whose reads face each other. Also without the flag were 81% of the pairs 20 kb or more apart and all inter-chromosomal pairs.
   - Consequence: do not select reads with `samtools view -f 2`, and pass `-useBadMate` (or `--use_bad_mates` in Bis-SNP) to tools that test this flag. Measure distances from POS and PNEXT.
2. The usage string calls the second argument `outputPrefix`; it is the BAM path.
3. `-clip5` and `-clip3` are BWA clipping penalties, not numbers of bases to trim, and their default 0 replaces BWA's default of 5. `-score 0` replaces BWA's default minimum score of 30.
4. `-h` does nothing. Usage is printed when the class is run with fewer than 4 tokens, and a usage error exits with status 0.
5. The pair choice is order-dependent in two places: the "both ends MAPQ above 0" test and the "same contig preferred" test favour only the later candidate, and a complete tie keeps the first pair seen. The `-enzymeList` step itself is symmetric.
6. `-enzymeList` acts only when MAPQ sums are equal, which is mostly at MAPQ 0. On one cell it changed the placement of 4.7% of pairs, 97% of them with an end below MAPQ 30, and lowered the pair alignment score in 84% of the changed pairs. It has almost no effect on reads that pass a MAPQ 30 filter.
7. The log line "Unique mapped and paired" counts every pair written; there is no uniqueness or MAPQ test behind it.
8. `-rgId` without `-rgSm` writes `RG` tags without an `@RG` header line.
9. `-buffer 0` divides by zero; a batch holds `buffer + 1` pairs.
10. `NM` and `MD` are in converted space until `samtools calmd` is run.
11. One alignment is kept per read end. The rest of a chimeric read stays soft-clipped and is not realigned. Measured on four scNOMe-HiC cells, realigning those parts would add 0.4% to 1.1% to the bases used for methylation calling.

`BhmemMultRef`, `BsmemMultRef` (untested classes)

12. `ambiguous.bam` receives the same records as `merged.bam`; there is no test for ambiguity.
13. The lift-over changes POS only. The contig name is set back to the original and the CIGAR is not adjusted.
14. Genome choice: complete ties favour the first `-refs` genome, and `-ratio` is a difference of MAPQ sums, not a ratio. With the default 0, a pair with equal MAPQ sums in both genomes goes to neither `alignToRef` file.
15. `-noHicMode`, `-noChimericMode` (`BhmemMultRef`) and `-chimericMode` (`BsmemMultRef`) are parsed and never read.
16. A `-refs` value without a directory component gives the index path `null/Bisulfite_Genome/...`.
17. Missing `-refs` or `-chains` ends in a stack trace instead of the usage text. `BsmemMultRef` prints `BhmemMultRef`'s usage string. No `XG`/`XR` tags are written.

Installation and wrapper

18. `jbwa.jar` and `libbwajni.so` are not in git; a fresh clone must run `install.sh` (or `scripts/quest_install.sh`) before `Bhmem` can start. The fat jar is in git.
19. `install.sh` links `bisulfitehic-default.jar` to the newest entry of `target/`, which can be a directory, and never corrects an existing link. `bisulfitehicMap` depends on that link and on `BISHIC`, and cannot pass `-pbat`, `-nonDirectional`, `-outputMateDiffChr` or `-enzymeList`. Call `java` directly instead.
20. `install.sudo.sh` sets `BISHIC` to `<repo>/target` when it was empty.
