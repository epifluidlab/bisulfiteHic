# bisulfiteHic usage guide

This guide is for someone who wants to use the package from another project without reading the source. It says what each tool takes, what it writes, and which options matter. The README has worked examples for the published analyses; this guide is the reference.

- [1. What is in the package](#1-what-is-in-the-package)
- [2. Install and run](#2-install-and-run)
- [3. Prepare a reference](#3-prepare-a-reference)
- [4. Standard workflow: FASTQ to methylation and contacts](#4-standard-workflow-fastq-to-methylation-and-contacts)
- [5. How to read a Bhmem BAM](#5-how-to-read-a-bhmem-bam)
- [6. Which tool for which question](#6-which-tool-for-which-question)
- [7. Conventions shared by the Java analysis tools](#7-conventions-shared-by-the-java-analysis-tools)
- [8. Tests](#8-tests)
- [9. Before you trust a result](#9-before-you-trust-a-result)

Reference pages, one per tool group:

| Page | Covers |
|---|---|
| [aligners.md](aligners.md) | `Bhmem`, `BhmemMultRef`, `BsmemMultRef`, installer scripts, `bisulfitehicMap` |
| [concordance.md](concordance.md) | `MethyCorAcrossHiccups`, `MethyCorAcrossHiccupsSC`, `MethyCorEachChrScBw`, `MethyCorGlobal`, `MethyCorSepAlleleSC` |
| [allele_specific.md](allele_specific.md) | `LongRangeAsm`, `AlleImbalanceExpr`, `SepAlleleMethySC`, `SepAlleleMethySC2`, `EnhanceHicMatrix` |
| [matrix_statistics.md](matrix_statistics.md) | `SingleBaseResHicMatrix`, `SingleBaseResHicMatrixV2`, `SingleCellTest`, `TestIndepPoiBinom`, `TestIndepPoiBinomOnBam` |
| [scripts.md](scripts.md) | Python, Perl and R helper scripts |

Each tool section has the same parts: purpose, command, inputs, options (every option, with its default and what the code does), output (column by column), behaviour notes, pitfalls and a minimal example. Each page ends with a list of known limitations.

## 1. What is in the package

The package handles bisulfite-converted Hi-C data: Methyl-HiC, NOMe-HiC, and their single-cell versions. One read pair carries both a chromatin contact and the methylation state of each end.

| Part | What it does |
|---|---|
| Aligner (`Bhmem`) | Aligns bisulfite Hi-C read pairs and writes one BAM. This is the part most projects need. |
| Java analysis tools (18 classes in one jar) | Read-level methylation concordance between two loci, allele-specific methylation, concordance matrices and tests. |
| Helper scripts | Demultiplexing, BAM to Juicer text, read-pair QC summaries, methylation track smoothing and tests, TAD calling. |

Methylation and SNP calling are not in this package. The workflow uses Bis-SNP for that (README, Analysis 2).

## 2. Install and run

**What a fresh clone contains.** The analysis jar `target/bisulfitehic-0.38-jar-with-dependencies.jar` is in git, and so are the library jars under `lib/`. The BWA bridge is not: `jbwa/jbwa-1.0.0/jbwa.jar` and `jbwa/jbwa-1.0.0/src/main/native/libbwajni.so` are built by `install.sh`. So the analysis tools run from a clone as it is, and `Bhmem` runs only after the install step.

**Build.** `install.sh` builds the BWA bridge with `make` and the jar with Maven. It needs `JAVA_HOME` set to a JDK 8 (the jbwa build uses `javah`), `gcc`, `make`, `wget`, zlib and Maven. On the Quest cluster, `scripts/quest_install.sh` does a full install with the right modules and runs the regression tests. Details: [aligners.md, Installing and running](aligners.md#installing-and-running).

**Java to run.** The classes are Java 8 bytecode. The tests run them on Java 11 and 17, and the sc_NOMeHiC pipeline runs `Bhmem` on Java 17.

**Running the aligner.** `<repo>` is the checkout directory.

```
java -Xmx20G \
  -Djava.library.path=<repo>/jbwa/jbwa-1.0.0/src/main/native \
  -cp "<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar:<repo>/jbwa/jbwa-1.0.0/jbwa.jar" \
  main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem \
  <ref.fa> <out.bam> <R1.fq.gz> <R2.fq.gz> [options]
```

**Running an analysis tool.** One classpath works for all of them (it is the one the unit tests use):

```
java -Xmx20G \
  -cp "<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar:<repo>/lib/*" \
  main.java.edu.mit.compbio.bisulfitehic.utils.<ClassName> <positional arguments> [options]
```

**Getting the option list.** Run the class with no arguments. `-h` does not work in any of the tools: it is accepted and ignored.

**Python scripts** need Python 3 with `pysam`, `numpy`, `scipy`, `bx-python`, `pyBigWig` and `tqdm`. See [scripts.md](scripts.md).

Do not use the `bisulfitehicMap` wrapper. It cannot pass the options that single-cell data needs, and it depends on a symlink that `install.sh` can create wrongly.

## 3. Prepare a reference

`Bhmem` needs these files next to the reference FASTA (here `genome.fa`):

| File | Made with |
|---|---|
| `Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa` and `Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa` | `bismark_genome_preparation <directory of genome.fa>` |
| BWA index of each of the two converted FASTA files | `bwa index <converted.fa>`, once for each |
| `genome.dict`, or `genome.fa.dict`, or `genome.fa.fai` | `samtools dict genome.fa -o genome.dict`, or `samtools faidx genome.fa` |

The file names under `Bisulfite_Genome` are fixed in the code. A `bwa` binary is built together with the BWA bridge under `<repo>/jbwa/jbwa-1.0.0/bwa-*/bwa`.

Optional, for `-enzymeList`: a tab-separated BED file of restriction sites (contig, start, end), plain or `.gz`. Lines starting with `#` are skipped; a blank line or a `track` line stops the run.

## 4. Standard workflow: FASTQ to methylation and contacts

These are the steps the sc_NOMeHiC pipeline runs for each cell. For bulk libraries, skip step 1 and change the `Bhmem` options as shown in step 3.

**Step 1. Demultiplex (single-cell plates only).** Splits a plate FASTQ pair by the 6-base barcode at the start of R1.

```
perl <repo>/src/perl/demultiplex_ecker_scWGBS.pl <out_prefix> <barcode_file> <R1.fastq.gz> <R2.fastq.gz>
```

**Step 2. Trim.** Not part of this package. The pipeline uses `trim_galore --paired_end` with `--clip_R1`, `--clip_R2`, `--three_prime_clip_R1` and `--three_prime_clip_R2`.

**Step 3. Align.** Choose the options by library type:

| Library | Options | In-silico conversion |
|---|---|---|
| Bulk Methyl-HiC or NOMe-HiC (directional) | `-outputMateDiffChr -buffer 100000 -t <threads> -rgId <id> -rgSm <sample>` | R1 C to T, R2 G to A |
| Single-cell Methyl-HiC or scNOMe-HiC (PBAT-type) | `-nonDirectional -pbat -outputMateDiffChr -buffer 1000 -enzymeList <sites.bed> -t <threads> -rgId <id> -rgSm <sample>` | R1 G to A, R2 C to T |

- Without `-outputMateDiffChr`, all contacts between chromosomes are discarded.
- Give `-rgId` and `-rgSm` together.
- The second positional argument is the BAM path, not a prefix.
- `-clip5` and `-clip3` are alignment penalties, not trimming lengths. Leave them alone.

**Step 4. Sort, fix mates, mark duplicates, recompute MD.** The raw BAM is unsorted, has no duplicate marks, and its `NM` and `MD` tags are in converted space.

```
samtools sort -n raw.bam | samtools fixmate -m - - | samtools sort - | samtools markdup - - \
  | samtools calmd -b - genome.fa > sample.calmd.bam
samtools index sample.calmd.bam
```

Every analysis tool in this package needs this BAM: coordinate-sorted, indexed, duplicate-marked, with `MD` tags from `calmd`.

**Step 5. QC summary.** Needs a name-sorted BAM.

```
samtools sort -n -o sample.byname.bam sample.calmd.bam
python <repo>/src/python/mh_reads_summary.v2.py --in_cram sample.byname.bam --out_summary sample.summary.txt
```

**Step 6. Methylation and SNPs.** Use Bis-SNP on `sample.calmd.bam` with `--use_bad_mates --mmq 30`, and `--nomeseq` for NOMe data (README, Analysis 2). The reads it uses are the non-duplicate reads with MAPQ of 30 or more.

**Step 7. Contacts.** Filter, convert to Juicer text, then build a `.hic` file.

```
samtools view -bh -q 30 -f 1 -F 1804 sample.calmd.bam > sample.good.bam
python <repo>/src/python/sam2juicer_new.py -s sample.good.bam -f <restriction_sites.txt> | gzip -c > sample.hic.txt.gz
java -jar juicer_tools.jar pre -d -q 30 -f <restriction_sites.txt> sample.hic.txt.gz sample.hic <genome>
```

`sam2juicer_new.py` does no filtering of its own.

## 5. How to read a Bhmem BAM

**What is in it.** One record for R1 and one for R2 of every read pair for which both ends aligned. There are no unmapped, secondary or supplementary records. Pairs with MAPQ 0 are written; filter downstream.

**One alignment per read end.** A read that crosses a ligation junction is written once, with the part that was not chosen soft-clipped. That part is not realigned.

**MAPQ** is BWA's value for the single-end alignment of the converted read to one converted genome.

**Strand flag (0x10).** The flag tells a methylation caller which bisulfite strand a read reports. Every tool in this package uses the same rule:

| Read | Flag 0x10 | Read reports |
|---|---|---|
| R1 | not set | C-to-T strand: C = methylated, T = unmethylated |
| R1 | set | G-to-A strand: G = methylated, A = unmethylated |
| R2 | set | C-to-T strand |
| R2 | not set | G-to-A strand |

With `-pbat` the aligner inverts the flag so that this rule also holds for PBAT-type libraries. SEQ and QUAL are always stored on the reference forward strand. The cost is that, in a `-pbat` BAM, the flag no longer gives the orientation in which the read aligned.

**Tags.** `XG` (`CT` or `GA`) is the converted genome the read aligned to. `XR` (`CT` or `GA`) is the in-silico conversion applied to the read. `AS` is the BWA score. `NM` and `MD` are in converted space until `samtools calmd` rewrites them. The analysis tools in this package read flags only, not `XG` or `XR`.

**TLEN and the proper-pair flag are not reliable. Do not filter on them.** In a `-pbat` BAM, pairs whose reads face each other look outward-facing by their flags. `samtools fixmate` then clears the proper-pair flag of those pairs. In one production cell, 55% of the usable pairs had no proper-pair flag, including most long-range contacts and all contacts between chromosomes.

- Do not use `samtools view -f 2`.
- Pass `-useBadMate` to the Java tools that have it, and `--use_bad_mates` to Bis-SNP.
- Measure the distance between the two ends from POS and PNEXT, not from TLEN.

Full details: [aligners.md](aligners.md#bhmem).

## 6. Which tool for which question

| Question | Tool | Main inputs | Page |
|---|---|---|---|
| Align bisulfite Hi-C reads | `Bhmem` | reference, R1, R2 | [aligners](aligners.md#bhmem) |
| Align to two parental genomes | `BhmemMultRef`, `BsmemMultRef` (untested) | two genomes, chain files | [aligners](aligners.md#bhmemmultref) |
| Split a single-cell plate by barcode | `demultiplex_ecker_scWGBS.pl` | plate FASTQ, barcode list | [scripts](scripts.md#perl-and-r-scripts) |
| How many reads are usable, cis, trans, long-range | `mh_reads_summary.v2.py` | name-sorted BAM | [scripts](scripts.md#hi-c-conversion-and-qc-scripts) |
| Make a `.hic` contact file | `sam2juicer_new.py` | filtered BAM, restriction sites | [scripts](scripts.md#hi-c-conversion-and-qc-scripts) |
| Call TADs from a `.hic` file | `call_topdom_on_hic.R` | `.hic` | [scripts](scripts.md#perl-and-r-scripts) |
| Is methylation at two loop anchors correlated on the same molecule | `MethyCorAcrossHiccups` | BAM, BEDPE | [concordance](concordance.md#methycoracrosshiccups) |
| The same across single cells | `MethyCorAcrossHiccupsSC`, `MethyCorEachChrScBw` | per-cell bigWigs | [concordance](concordance.md#methycoracrosshiccupssc) |
| Window-by-window concordance along a chromosome | `MethyCorGlobal` | BAM | [concordance](concordance.md#methycorglobal) |
| Allele-specific methylation linked over long range | `LongRangeAsm` | phased VCF, BAM | [allele-specific](allele_specific.md#longrangeasm) |
| Per-CpG methylation by allele in one cell | `SepAlleleMethySC`, `SepAlleleMethySC2` | VCF, BAM | [allele-specific](allele_specific.md#sepallelemethysc) |
| Allele-specific expression per gene | `AlleImbalanceExpr` | gene BED, VCF, RNA BAM | [allele-specific](allele_specific.md#alleimbalanceexpr) |
| CpG-by-CpG concordance table from read pairs | `SingleBaseResHicMatrix` | CpG list, BAM | [matrix and statistics](matrix_statistics.md#singlebasereshicmatrix) |
| Test independence of methylation between sites across cells | `SingleCellTest`, `TestIndepPoiBinom`, `TestIndepPoiBinomOnBam` | cell-by-site table or per-cell BAMs | [matrix and statistics](matrix_statistics.md#singlecelltest) |
| Smooth or test GCH/HCG methylation tracks | `smooth_mgch.py`, `binom_test_mgch.py`, `betabinom_test_mgch.py` | methylation BED | [scripts](scripts.md#gchhcg-methylation-track-scripts) |

## 7. Conventions shared by the Java analysis tools

**Input BAM.** Coordinate-sorted and indexed, duplicates flagged, `MD` tag present (step 4 above). The tools rebuild the reference bases from `MD`; they never open a reference FASTA. Two tools need a name-sorted BAM instead and say so in their sections (`SepAlleleMethySC2`, `TestIndepPoiBinomOnBam`).

**Default read filters.** Unmapped, secondary, duplicate and QC-failed records are dropped. Minimum MAPQ is 30 or 20 depending on the tool (`-minMapQ`). Minimum base quality for a methylation call is 5 in most tools (`-minBaseQ`). Records without the proper-pair flag are dropped unless `-useBadMate` is given, and for Hi-C data it should be given (section 5). Supplementary records are not filtered, which does not matter for `Bhmem` BAMs because they contain none.

**Methylation context.** Tools with pattern options take an IUPAC pattern and the 1-based position of the cytosine in it. The pattern is matched against the reference, not the read.

| Data | Call | Options |
|---|---|---|
| Methyl-HiC, CpG methylation | CG | defaults: `-methyPatternSearch CG -methyPatternCytPos 1 -bsConvPatternSearch CH -bsConvPatternCytPos 1` |
| NOMe-HiC, accessibility | GCH | `-methyPatternSearch GCH -methyPatternCytPos 2 -bsConvPatternSearch WCH -bsConvPatternCytPos 2` |
| NOMe-HiC, endogenous methylation | HCG or WCG | `-methyPatternSearch HCG -methyPatternCytPos 2 -bsConvPatternSearch WCH -bsConvPatternCytPos 2` |

Tools without pattern options are fixed to CpG with a `CH` conversion filter. On NOMe data that filter counts methylated GCH sites as unconverted, which is why the NOMe examples use `WCH`. These tools cannot be switched to `WCH`, so check them before using them on NOMe data: `MethyCorGlobal`, `MethyCorSepAlleleSC`, `SepAlleleMethySC`, `SepAlleleMethySC2`, `SingleBaseResHicMatrix`, `SingleBaseResHicMatrixV2`, `TestIndepPoiBinomOnBam`.

**Bisulfite conversion filter.** On by default. A read is dropped when 40% or more of its conversion-context sites are unconverted. `-turnOffBisulfiteFilter` switches it off where the option exists.

**Distance between the two ends** is `|POS - PNEXT|` on the same contig. TLEN is not used.

**Alleles.** In the allele-specific tools "ref" and "alt", or A and B, are the first and second allele of the genotype in the VCF, not REF and ALT.

**Output.** Most tools write gzip-compressed, tab-separated text without a header line. Missing values are `NaN` or `NA`, depending on the tool.

**Command line.** Positional arguments and options can be mixed. List options are repeated: `-bams a.bam -bams b.bam`. A usage error prints the usage text and exits with status 0, so check that the output file exists.

## 8. Tests

```
bash test/regression/run_tests.sh   # Bhmem on a synthetic reference, about one minute
bash test/unit/run_all.sh           # analysis classes and helper scripts on synthetic inputs
```

Not covered by any test: `BhmemMultRef`, `BsmemMultRef`, `MethyCorAcrossHiccupsSC`, `MethyCorSepAlleleSC`, `mh_reads_summary.v2.py`, `mh_hic_filter.py`, `eigen_adjust_direction.py`, `sam2juicer.py`, `call_topdom_on_hic.R`.

## 9. Before you trust a result

The path the sc_NOMeHiC pipeline uses (`Bhmem`, the samtools steps, `mh_reads_summary.v2.py`, `sam2juicer_new.py`, `demultiplex_ecker_scWGBS.pl`) is tested, and its bugs found in November 2026 are fixed (README, change log). The other tools were written for specific published analyses. Many have options that are parsed and ignored, and some have behaviour that looks unintended. These are listed at the end of each reference page. The ones most likely to change a result:

1. **Proper-pair flag and TLEN in `-pbat` BAMs** (section 5). Tools that require the flag lose most Hi-C contacts. `SingleBaseResHicMatrix` and its V2 always require it. ([aligners](aligners.md#known-limitations))
2. **CpG list of `SingleBaseResHicMatrix`:** the start column must be the 1-based position of the C. A standard BED gives an empty result. ([matrix and statistics](matrix_statistics.md#known-limitations))
3. **`-randomSampleReads` in `MethyCorAcrossHiccups`** repeats the shuffle; it does not subsample. Each pair is written up to that many times. ([concordance](concordance.md#known-limitations))
4. **`TestIndepPoiBinomOnBam -chr`** only labels the output; reads of all chromosomes are pooled. ([matrix and statistics](matrix_statistics.md#known-limitations))
5. **Insertions and deletions** shift the coordinates of later methylation calls in the allele-specific and single-CpG tools. ([allele-specific](allele_specific.md#known-limitations))
6. **`BhmemMultRef` and `BsmemMultRef`:** `ambiguous.bam` is a copy of `merged.bam`, and the lift-over changes POS only. ([aligners](aligners.md#known-limitations))
7. **Permutation and shuffle modes are unseeded** (except `MethyCorEachChrScBw -seed`), so their output differs between runs.
