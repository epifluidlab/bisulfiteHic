# Allele-specific tools

Part of the [usage guide](USAGE.md). These five Java tools split reads by the allele they carry at heterozygous SNPs.

| Tool | Input | Output | Choose it when |
|---|---|---|---|
| [`LongRangeAsm`](#longrangeasm) | phased VCF + BAM | per-SNP methylation of each allele: at the SNP read, at its distant mate, and locally | you want allele-specific methylation linked over long range (README Analysis 5) |
| [`AlleImbalanceExpr`](#alleimbalanceexpr) | gene BED + VCF + BAM | per-gene allele read counts and test | allele-specific expression from an RNA BAM |
| [`SepAlleleMethySC`](#sepallelemethysc) | VCF + coordinate-sorted BAM | per-CpG methylation by allele | one cell, reads overlapping a SNP |
| [`SepAlleleMethySC2`](#sepallelemethysc2) | VCF + name-sorted BAM | per-CpG methylation by allele | as above, also assigning the mate of a SNP read |
| [`EnhanceHicMatrix`](#enhancehicmatrix) | CpG list + BAM (+ BEDPE) | methylation-concordance statistics per window pair | window-level concordance with a permutation test |

- The unit tests (`test/unit/allele/`) use C-to-T strand reads only.
- Two helpers live in `lib/ccInference-0.19-slim.jar` without source in this repository (the conversion filter and the IUPAC matcher); what is said about them was read from the compiled classes.
- Sentences that start with "Not verified:" rest on library behaviour that was not checked.

---

## Shared behaviour (applies to all five tools unless a section says otherwise)

**Classpath.** The README runs these tools with three jars, in this order:

```
java -Xmx20G -cp "/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.<ClassName> ...
```

The fat jar contains the tool classes, htsjdk, args4j, commons-math3 and log4j, but it does not contain `CcInferenceUtils` (ccInference jar) or `GenomeLocParser` (GATK jar), so all three jars are needed. `-Djava.library.path` (jbwa) is not used by these tools in the README.

**Options and positional arguments.** Options are args4j options with a single dash. The README and the unit tests put the options after the positional arguments. Each tool checks `args.length` (all tokens, options included) against the number of positional arguments before parsing; if there are fewer tokens it prints the usage string and the option list to stderr and returns without a non-zero exit call.

**`-h`.** In all five tools the `help` field is tested before the arguments are parsed, so `-h` is never seen as true. The option list is printed only when too few tokens are given (for example, when the class is run with no arguments, which is what the README recommends).

**Strand classes used below.** "C->T strand read" means first-of-pair on the forward strand or second-of-pair on the reverse strand. "G->A strand read" means first-of-pair on the reverse strand or second-of-pair on the forward strand. A read without the paired flag is treated as first-of-pair.

**How a methylation call is made (pattern options).** The reference bases under a read are rebuilt from the read itself with htsjdk `SequenceUtil.makeReferenceFromAlignment`, which needs the `MD` tag (it throws `SAMException: Cannot create reference from SAMRecord with no MD tag` otherwise). No reference FASTA is read. Soft-clipped bases are removed. For a G->A strand read the read, the reference and the qualities are reverse-complemented first, so the cytosine is always looked for as C/T.
- `-methyPatternSearch` is a pattern matched against the reference bases, with IUPAC codes allowed in the pattern (standard IUPAC sets, e.g. H = A/C/T, W = A/T). A reference base that is not A/C/G/T never matches.
- `-methyPatternCytPos` is the 1-based position inside the pattern whose read base is inspected. Read base C with base quality >= `-minBaseQ` is a methylated call, read base T with base quality >= `-minBaseQ` is an unmethylated call, anything else gives no call. The code does not check that the pattern letter at that position is a C.
- `-methyPattern2ndSearch` / `-methyPattern2ndCytPos` (LongRangeAsm only) run the same search a second time with another pattern and fill extra columns.
- `SepAlleleMethySC` and `SepAlleleMethySC2` have no pattern options; they use a fixed reference `CG` test.

**Bisulfite conversion filter (`-bsConvPatternSearch`, `-bsConvPatternCytPos`, `-turnOffBisulfiteFilter`).** Unless `-turnOffBisulfiteFilter` is given, each read is passed to `CcInferenceUtils.bisulfiteIncompleteReads(read, ref, negStrand, secondPair, false, true, true, 0.1, 0.4, 0, pattern, cytPos)`. With these fixed arguments the method counts (a) the windows of the read, oriented to the cytosine strand, that literally match the conversion pattern (an unconverted cytosine in that context) and (b) the windows of the reference that match the pattern. The read is rejected when (a)/(b) >= 0.4, or when (b) is 0 and (a) is not 0. A rejected read gives no methylation calls. With the fixed arguments the 0.1 value and the cytosine position are not used for the decision, so `-bsConvPatternCytPos` does not change the result. `SepAlleleMethySC` and `SepAlleleMethySC2` always use pattern `CH`.

**Coordinates of calls.** Calls are indexed by their offset in the soft-clip-trimmed read and converted to a coordinate as `alignmentStart + offset`. Deleted reference bases have no entry in the arrays and inserted bases do, so in a read with an insertion or deletion the coordinates of calls after the indel are shifted by the indel length.

**Logging.** Each tool calls log4j `BasicConfigurator.configure()` and logs progress and the running time.

---

## `LongRangeAsm`

**Purpose** — For every heterozygous SNP it collects the read pairs whose two ends are far apart (default at least 20 kb), assigns each pair to allele A or B from the base at the SNP, and reports per-SNP counts of methylated and total calls on the SNP-side read and on the distal mate, for each allele. Choose it when you want one aggregated row per SNP with a distance filter and a configurable methylation pattern (CG, GCH, HCG ...). `SepAlleleMethySC`/`SepAlleleMethySC2` instead write one line per call with a fixed CG pattern and no distance filter.

**Command**

```
java -Xmx20G -cp "<fat jar>:<ccInference-0.19-slim.jar>:<gatk-package-distribution-3.3.jar>" \
  main.java.edu.mit.compbio.bisulfitehic.utils.LongRangeAsm [opts] genotype.vcf.gz input.bam output.ASM.tab.gz output.ref.tab.gz output.alt.tab.gz
```

Positional arguments (at least 5 tokens are required):
1. `genotype.vcf.gz` — SNP VCF.
2. `input.bam` — alignments.
3. `output.ASM.tab.gz` — main per-SNP table.
4. `output.ref.tab.gz` — per-cytosine table of allele A reads.
5. `output.alt.tab.gz` — per-cytosine table of allele B reads.

**Inputs**
- VCF: the tabix index (`<vcf>.tbi`) must exist and be readable, otherwise the tool throws `UnsortedFileException`. The sequence dictionary is taken from the VCF itself. Records that are filtered, indels, mixed or not point events are skipped.
  - Default mode: the genotype of sample number `-sample` (0-based) must be not filtered, heterozygous and phased. Allele A is the first allele of the phased genotype and allele B the second (for `1|0` allele A is the ALT base). A second record at a position already loaded is ignored.
  - `-partialVcf`: genotypes are not read. Allele A is the REF base and allele B is the first ALT base. `-sample` is not used.
- BAM: the file name must end in `.bam`. It is read with region queries and mate queries, so it must be coordinate-sorted and indexed. Reads need the `MD` tag (the README pipeline adds it with `samtools calmd`). The mate fields (mate contig and mate start) must be filled, because the end-to-end distance is computed from them. The mate record itself is fetched with `queryMate`.
- Chromosome names: the SNP contig name from the VCF is used unchanged to query the BAM, so both files must use the same names. The README example uses `-region 22` (no `chr` prefix).
- Pattern options: see "Shared behaviour".

**Options**

| option | default | what the code does |
|---|---|---|
| `-region` | null | Restricts the SNPs that are loaded. `chr:start-end` queries that VCF interval; `chr` alone queries the whole contig (the contig must be in the VCF dictionary). Null loads all SNPs. It does not restrict where mates may map. |
| `-turnOffBisulfiteFilter` | false | Skips the bisulfite conversion filter for both ends. |
| `-useBadMate` | false | Reads without the proper-pair flag are accepted (in the tool's own filter and in the methylation helper). |
| `-diffChr` | false | Pairs whose mate is on another contig are accepted (they are also exempt from the proper-pair test). The per-cytosine tables are then written once at the end of the run. |
| `-useUnpaired` | false | Removes the "read must be paired" test in the tool's own filter only. The methylation helper still rejects unpaired reads, so no unpaired read is ever counted. See Pitfalls. |
| `-methyPatternSearch` | `CG` | Primary reference pattern for methylation calls. |
| `-methyPatternCytPos` | 1 | 1-based position of the inspected base in `-methyPatternSearch`. |
| `-methyPattern2ndSearch` | null | Second pattern. When set, 18 extra columns are written. The usage text says "default: CG" but the default is null (off). |
| `-methyPattern2ndCytPos` | 1 | 1-based position of the inspected base in `-methyPattern2ndSearch`. |
| `-bsConvPatternSearch` | `CH` | Pattern used by the bisulfite conversion filter. |
| `-bsConvPatternCytPos` | 1 | Passed to the filter, but it does not change the filter result (see "Shared behaviour"). |
| `-minMapQ` | 30 | Minimum mapping quality of the SNP read and of its mate. |
| `-minBaseQ` | 20 | Minimum base quality of the base at the SNP and of every base used for a methylation call. (The usage text wrongly says "mapping quality".) |
| `-minDist` | 20000 | Minimum distance between the two ends on the same contig. |
| `-maxDist` | 1000000000 | Maximum distance between the two ends on the same contig. |
| `-coverageRef` | 5 | Minimum `AlleleA_local_count` for a SNP to be written. It applies to allele A, which is not necessarily the REF allele. |
| `-coverageAlt` | 5 | Minimum `AlleleB_local_count` for a SNP to be written. |
| `-qual` | 30 | Parsed but never used. No genotype-quality filter exists in the code. |
| `-partialVcf` | false | Use REF/ALT of the record instead of the phased genotype (see Inputs). |
| `-sample` | 0 | 0-based index of the genotype column used in default mode. |
| `-h` | false | No effect (see "Shared behaviour"). |

**Output**

1. Main table (argument 3). Gzip-compressed, tab-separated, one header line starting with `#chr`, then one row per SNP that passes the coverage thresholds, in the order of the sorted SNP map (`TreeMap` keyed by GATK `GenomeLoc`). `*_C` = number of methylated calls, `*_CT` = number of all calls (methylated + unmethylated), `*_meth_mean` = `*_C / *_CT`, printed as `NaN` when `*_CT` is 0.

| col | header label | meaning |
|---|---|---|
| 1 | `#chr` | SNP contig |
| 2 | `start` | SNP position (1-based, as in the VCF) |
| 3 | `end` | SNP end (equal to column 2 for a SNP) |
| 4 | `AlleleA` | allele A base |
| 5 | `AlleleA_count` | allele A reads at the SNP whose mate record was found and passes the read filter |
| 6 | `AlleleA_C` | methylated calls in the SNP-side reads counted in column 5 |
| 7 | `AlleleA_CT` | all calls in those reads |
| 8 | `AlleleA_meth_mean` | column 6 / column 7 |
| 9 | `AlleleB` | allele B base |
| 10 | `AlleleB_count` | as column 5 for allele B |
| 11 | `AlleleB_C` | as column 6 for allele B |
| 12 | `AlleleB_CT` | as column 7 for allele B |
| 13 | `AlleleB_meth_mean` | column 11 / column 12 |
| 14 | `End2_anchor_point` | alignment starts of the mates, joined with `:`; a mate start within 500 bp of one already listed is not added; no contig name; unordered; empty when no mate was usable |
| 15 | `End2_AlleleA_C` | methylated calls in the mates of allele A reads |
| 16 | `End2_AlleleA_CT` | all calls in the mates of allele A reads |
| 17 | `End2_AlleleA_meth_mean` | column 15 / column 16 |
| 18 | `End2_AlleleB_C` | methylated calls in the mates of allele B reads |
| 19 | `End2_AlleleB_CT` | all calls in the mates of allele B reads |
| 20 | `End2_AlleleB_meth_mean` | column 18 / column 19 |
| 21 | `MaxEnd2Dist` | largest `abs(mate start - SNP position)` over the anchors in column 14; -1 when there is none |
| 22 | `AlleleA_local_count` | allele A reads at the SNP, whether or not the mate record was found/passed |
| 23 | `AlleleA_local_C` | methylated calls in the reads counted in column 22 |
| 24 | `AlleleA_local_CT` | all calls in those reads |
| 25 | `AlleleA_local_meth_mean` | column 23 / column 24 |
| 26 | `AlleleB_local_count` | as column 22 for allele B |
| 27 | `AlleleB_local_C` | as column 23 for allele B |
| 28 | `AlleleB_local_CT` | as column 24 for allele B |
| 29 | `AlleleB_local_meth_mean` | column 27 / column 28 |

Only with `-methyPattern2ndSearch` (same definitions, second pattern):

| col | header label | meaning |
|---|---|---|
| 30 | `AlleleA_2nd_local_C` | 2nd pattern, allele A, SNP-side reads, mate not required |
| 31 | `AlleleA_2nd_local_CT` | |
| 32 | `AlleleA_2nd_local_meth_mean` | |
| 33 | `AlleleB_2nd_local_C` | 2nd pattern, allele B, SNP-side reads, mate not required |
| 34 | `AlleleB_2nd_local_CT` | |
| 35 | `AlleleB_2nd_local_meth_mean` | |
| 36 | `AlleleA_2nd_C` | 2nd pattern, allele A, SNP-side reads with a usable mate |
| 37 | `AlleleA_2nd_CT` | |
| 38 | `AlleleA_2nd_meth_mean` | |
| 39 | `AlleleB_2nd_C` | 2nd pattern, allele B, SNP-side reads with a usable mate |
| 40 | `AlleleB_2nd_CT` | |
| 41 | `AlleleB_2nd_meth_mean` | |
| 42 | `End2_AlleleA_2nd_C` | 2nd pattern, mates of allele A reads |
| 43 | `End2_AlleleA_2nd_CT` | |
| 44 | `End2_AlleleA_2nd_meth_mean` | |
| 45 | `End2_AlleleB_2nd_C` | 2nd pattern, mates of allele B reads |
| 46 | `End2_AlleleB_2nd_CT` | |
| 47 | `End2_AlleleB_2nd_meth_mean` | |

The table has 29 columns without and 47 columns with `-methyPattern2ndSearch`.

Check of the README "Analysis 5" R script: it reads `x[,c(8,13,17,20,25,29)]` as the six primary-pattern means, and the count pairs `(6,7)`, `(11,12)`, `(15,16)`, `(18,19)`, `(23,24)`, `(27,28)`, `(30,31)`, `(33,34)`, `(36,37)`, `(39,40)`, `(42,43)`, `(45,46)` as (methylated, total) for allele A, allele B, end-2 A, end-2 B, local A, local B and the six second-pattern groups in the order of columns 30-47. These indices are consistent with the code, and the script's use of `s[2]-s[1]` as the unmethylated count is consistent with `*_CT` being a total. The script needs the 47-column layout, which the README command produces because it passes `-methyPattern2ndSearch HCG`.

2. Per-cytosine tables (arguments 4 and 5; argument 4 = allele A reads, argument 5 = allele B reads). Gzip-compressed, tab-separated, no header, rows in no particular order, primary pattern only.

| col | meaning |
|---|---|
| 1 | contig of the read that carries the call (the mate's own contig for a trans mate) |
| 2 | position - 1 |
| 3 | position (1-based). For a C->T strand read this is the cytosine. For a G->A strand read it is the forward-strand base paired with the cytosine (for `CG`: the G). |
| 4 | column 5 / column 6 |
| 5 | number of methylated calls at this position |
| 6 | number of all calls at this position |

**Behaviour notes**
- Reads fetched at a SNP are dropped when they are unmapped, not primary, below `-minMapQ`, vendor-failed, duplicate, unpaired, or not proper pair (unless `-useBadMate`, or `-diffChr` with the mate on another contig).
- Distance: `abs(alignment start - mate alignment start)` from the read's own mate fields, never TLEN. On the same contig the read is kept when `minDist <= distance <= maxDist`. A mate on another contig or an unmapped mate is dropped unless `-diffChr`.
- One record per read name is used per SNP (the first one returned by the query).
- The mate record is fetched with `queryMate` and must pass the same filter. If it is missing or fails, the SNP-side read still counts in the `*_local_*` columns only.
- Allele assignment: the base at the SNP is located through the CIGAR (soft clips handled); the read is skipped when it has no base there or the base quality is below `-minBaseQ`. For a C/T SNP (either order) C->T strand reads are skipped. For a G/A SNP G->A strand reads are skipped. Otherwise, on a C->T strand read a read base T also matches an allele C, and on a G->A strand read a read base A also matches an allele G; all other comparisons are exact. Allele A is tested first; a read matching neither allele adds nothing to the allele counts.
- A SNP-side read with no primary-pattern call (none in the read, or rejected by the bisulfite filter) is skipped completely: it is not counted in any `*_count` column, even when it has second-pattern calls.
- A SNP row is written when `AlleleA_local_count >= coverageRef` and `AlleleB_local_count >= coverageAlt` (defaults 5 and 5) and at least one local `*_CT` is non-zero. There is no threshold on the end-2 columns.
- The "local" columns are not short-range reads: every read in them already passed the distance filter through its own mate fields.

**Pitfalls**
- By default a pair must carry the proper-pair flag and be at least 20 kb apart. The README example passes `-useBadMate`.
- `-qual` and `-h` do nothing. `-bsConvPatternCytPos` does not change the result.
- `-useUnpaired` cannot add reads. With `-useUnpaired` alone, an unpaired read that reaches the proper-pair test makes htsjdk throw `IllegalStateException: Inappropriate call if not paired read` (htsjdk `SAMRecord.getProperPairFlag`); with `-useBadMate` it fails the distance test; with `-diffChr` it passes the tool's filter and is then rejected by the helper.
- `-coverageRef`/`-coverageAlt` and the "ref"/"alt" output files refer to allele A / allele B. Without `-partialVcf` allele A is the first phased allele, which is the ALT base for a `1|0` genotype.
- The per-cytosine tables are filled before the coverage test, so they also contain reads of SNPs that are not in the main table and reads without a usable mate. A read (and its mate) that overlaps several SNPs is added once per SNP.
- A methylation call at the SNP position itself is not excluded.
- With `-diffChr`, `End2_anchor_point` and `MaxEnd2Dist` mix coordinates of different contigs (the mate start is compared with the SNP position without a contig check).
- With `-partialVcf`, two records at the same position stop the run (`GenomeLocSortedSet already contains the GenomeLoc`, the GATK jar); the default mode skips the second record.
- If fewer than five positional arguments are given but at least five tokens in total, the tool fails with an index error instead of printing the usage.
- Reads with indels: see "Coordinates of calls".

**Minimal example** (README, Analysis 5)

```
java -Xmx20G -cp "/your/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/your/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.LongRangeAsm /your/path/to/phased_vcf_files/testOut.vcf.gz testOut.b37.calmd.bam testOut.ASM_GCH_HCG_20kb_minCov2.chr22.tab.gz testOut.ASM_GCH_HCG_20kb_minCov2.ref.chr22.tab.gz testOut.ASM_GCH_HCG_20kb_minCov2.alt.chr22.tab.gz -minDist 20000 -sample 0 -coverageRef 2 -coverageAlt 2 -methyPatternSearch GCH -methyPatternCytPos 2 -methyPattern2ndSearch HCG -methyPattern2ndCytPos 2 -bsConvPatternSearch WCH -bsConvPatternCytPos 2 -useBadMate -minMapQ 30 -region 22
```

---

## `AlleImbalanceExpr`

**Purpose** — For every heterozygous SNP inside a gene interval it counts the reads carrying allele A, allele B, the REF base or another base, and gives a two-sided binomial p-value for A versus B (p = 0.5), per SNP or summed per gene. It compares read bases with the alleles exactly (no bisulfite handling) and makes no methylation calls, so it is the tool for allelic imbalance in ordinary (for example RNA-seq) alignments; the other tools in this part are for bisulfite reads.

**Command**

```
java -Xmx20G -cp "<fat jar>:<ccInference-0.19-slim.jar>:<gatk-package-distribution-3.3.jar>" \
  main.java.edu.mit.compbio.bisulfitehic.utils.AlleImbalanceExpr [opts] gene.bed[.gz] genotype.vcf.gz input.bam output.ASE.tab.gz
```

Positional arguments (at least 4 tokens are required):
1. `gene.bed[.gz]` — gene intervals.
2. `genotype.vcf.gz` — SNP VCF.
3. `input.bam` — alignments.
4. `output.ASE.tab.gz` — output table.

**Inputs**
- Gene BED: tab-separated; read through gzip when the name ends in `.gz`; lines starting with `#` are skipped. Columns 1-3 (contig, 0-based start, end) define the interval; the SNPs tested are those at 1-based positions `start+1 .. end`. The whole input line is copied to the output. The output header assumes six BED columns (`chr, start, end, name, ., strand`), but the code does not check the column count. The contig names are turned into GATK `GenomeLoc` objects with the VCF sequence dictionary. Two genes with identical coordinates are both reported.
- VCF: same rules as `LongRangeAsm` (tabix index required; filtered, indel, mixed and non-point records skipped; default mode needs a not-filtered, heterozygous, phased genotype in sample `-sample`; allele A = first phased allele, allele B = second; `-partialVcf` uses REF as A and the first ALT as B). The REF base is kept as a third value for the `AlleleRef` columns.
- BAM: file name must end in `.bam`; coordinate-sorted and indexed (region queries). No `MD` tag and no pairing information is needed. With `-waspFilter` the `vW` integer tag written by STAR is read.
- Chromosome names: BED and VCF must use the same names. The BAM is queried with the VCF name, or with `chr` + name when `-consistentPrefix` is given and the name does not start with `chr`.
- Methylation pattern options: none.

**Options**

| option | default | what the code does |
|---|---|---|
| `-chrom` | null | Keeps only BED rows whose contig equals this name (case-insensitive) and loads only this contig from the VCF (the name must be in the VCF dictionary). |
| `-mergedInSameGene` | false | Sums the counts of all SNPs of a gene and writes one row per gene with a different column layout (see Output). |
| `-minMapQ` | 30 | Minimum mapping quality of a read. |
| `-minBaseQ` | 5 | Minimum base quality of the base at the SNP. (The usage text wrongly says "mapping quality".) |
| `-minCov` | 5 | A row is written when the allele A count or the allele B count is at least this value (it is skipped only when both are below). |
| `-qual` | 30 | Parsed but never used. |
| `-partialVcf` | false | Use REF/ALT of the record instead of the phased genotype. |
| `-consistentPrefix` | false | Prepends `chr` to the contig name for the BAM query only, when the name does not already start with `chr`. It does not change BED or VCF names. |
| `-waspFilter` | false | Drops reads that have a `vW` tag whose value is not 1. Reads without the tag are kept. |
| `-sample` | 0 | 0-based index of the genotype column used in default mode. |
| `-h` | false | No effect (see "Shared behaviour"). |

**Output** — one gzip-compressed, tab-separated file (argument 4) with one header line starting with `#chr`. The first columns are the unchanged BED line; the numbers below assume a 6-column BED.

Default mode (one row per gene x SNP), 15 columns:

| col | header label | meaning |
|---|---|---|
| 1-6 | `#chr`, `gene_start`, `gene_end`, `gene_name`, `.`, `strand` | the BED line |
| 7 | `snp_pos` | SNP position (1-based) |
| 8 | `AlleleRef` | REF base of the VCF record |
| 9 | `AlleleRef_count` | reads whose base equals the REF base |
| 10 | `AlleleA` | allele A base |
| 11 | `AlleleA_count` | reads whose base equals allele A |
| 12 | `AlleleB` | allele B base |
| 13 | `AlleleB_count` | reads whose base equals allele B (and not allele A) |
| 14 | `Binomial_p` | two-sided binomial test of column 11 out of (column 11 + column 13) against 0.5, printed with 3 significant digits (`%.3g`) |
| 15 | `AlleleN_count` | reads whose base is neither allele |

`-mergedInSameGene` (one row per gene), 11 columns:

| col | header label | meaning |
|---|---|---|
| 1-6 | as above | the BED line |
| 7 | `snp_pos` | positions of all SNPs in the gene, joined with `:` |
| 8 | `AlleleA_count` | sum over the SNPs |
| 9 | `AlleleB_count` | sum over the SNPs |
| 10 | `Binomial_p` | binomial test on the summed counts |
| 11 | `AlleleN_count` | sum over the SNPs |

**Behaviour notes**
- Reads are dropped when unmapped, not primary, below `-minMapQ` (default 30), vendor-failed, or (with `-waspFilter`) tagged `vW` with a value other than 1. Duplicate-flagged reads are not dropped. Pairing, proper-pair flag and distance between ends are not examined.
- One record per read name is counted per SNP (the last one returned by the query), so two overlapping ends count once.
- The base at the SNP is located through the CIGAR (soft clips handled); a read without a base there or with base quality below `-minBaseQ` (default 5) is ignored.
- Allele assignment is an exact base comparison: allele A first, then allele B, otherwise the read goes to `AlleleN_count`. The REF count is independent of this.
- A SNP that lies in several genes is reported once per gene.

**Pitfalls**
- No bisulfite-aware matching: on bisulfite reads converted bases would be counted as the wrong allele or as `AlleleN`.
- `-qual` and `-h` do nothing.
- In `-mergedInSameGene` mode a read that covers two SNPs of the gene is counted once per SNP, and the summed A/B counts are only meaningful when allele A is the same haplotype at every SNP (phased genotypes; not guaranteed with `-partialVcf`).
- The header fits only a 6-column BED. A blank line in the BED file causes an index error.
- `-consistentPrefix` does not help when BED and VCF names differ from each other.
- The class comment at the top of the source describes an older column layout (see Known limitations).

**Minimal example** (shape from the usage string and `test/unit/allele/AlleImbalanceExprTest.java`; the README has no example)

```
java -Xmx20G -cp "/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.AlleImbalanceExpr genes.bed phased.vcf.gz sample.rna.bam sample.ASE.tab.gz
```

---

## `SepAlleleMethySC`

**Purpose** — For every heterozygous SNP it takes each read that overlaps the SNP, assigns it to allele A or B, and writes every CpG call of that read to an allele A or an allele B file, one line per call. Choose it for a coordinate-sorted, indexed BAM when only the read that covers the SNP should be used. `SepAlleleMethySC2` also assigns the mate; `LongRangeAsm` aggregates per SNP and filters by distance.

**Command**

```
java -Xmx20G -cp "<fat jar>:<ccInference-0.19-slim.jar>:<gatk-package-distribution-3.3.jar>" \
  main.java.edu.mit.compbio.bisulfitehic.utils.SepAlleleMethySC [opts] outputPrefix snp_vcf.gz input.bam
```

Positional arguments (at least 3 tokens are required):
1. `outputPrefix` — prefix of the two output files.
2. `snp_vcf.gz` — SNP VCF.
3. `input.bam` — alignments.

**Inputs**
- VCF: the tabix index must exist. Records that are filtered, indels, mixed or not point events are skipped. The genotype of the first sample (index 0, not configurable) must be not filtered and heterozygous. A phased genotype is not required by the code. Allele A is the first allele of the genotype as htsjdk returns it and allele B the second; they correspond to haplotypes only when the VCF is phased.
- BAM: coordinate-sorted and indexed (one region query per SNP). Reads need the `MD` tag. Reads must be paired.
- Chromosome names: the VCF contig name is used unchanged to query the BAM and is written to the output.
- Methylation pattern: fixed. A call is made at every reference `CG`; the bisulfite filter uses `CH`. There are no pattern options.

**Options**

| option | default | what the code does |
|---|---|---|
| `-minBaseQ` | 5 | Minimum base quality of the base at the SNP and of every base used for a call. |
| `-useBadMate` | false | Accepts reads without the proper-pair flag. |
| `-minMapQ` | 30 | Minimum mapping quality. |
| `-turnOffBisulfiteFilter` | false | Skips the bisulfite conversion filter. |
| `-h` | false | No effect (see "Shared behaviour"). |

**Output** — two gzip-compressed, tab-separated files without header: `<outputPrefix>.alleleA_methy.bedgraph.gz` and `<outputPrefix>.alleleB_methy.bedgraph.gz`. One line per call per read (calls are not merged per position; lines are not sorted).

| col | meaning |
|---|---|
| 1 | contig (from the VCF) |
| 2 | position - 1 |
| 3 | position (1-based) of the cytosine of the CpG on the forward strand, for reads of either strand |
| 4 | 1 = methylated, 0 = unmethylated |

**Behaviour notes**
- Reads are dropped when unmapped, not primary, below `-minMapQ` (30), vendor-failed, duplicate, unpaired, or not proper pair (unless `-useBadMate`). The distance between ends is not examined. The mate is not used.
- Allele assignment is the same as in `LongRangeAsm`: for a C/T SNP C->T strand reads are skipped, for a G/A SNP G->A strand reads are skipped, otherwise T matches an allele C on C->T strand reads and A matches an allele G on G->A strand reads; allele A is tested first. A read with no base at the SNP or with base quality below `-minBaseQ` (5) is skipped.
- A call: C->T strand read, reference `C` followed by `G`, read C = 1 and read T = 0. G->A strand read, reference `G` preceded by `C`, read G = 1 and read A = 0, reported at the position of the C.
- The call at the SNP position itself is not written.
- A read is processed once for every SNP it overlaps; both ends of a pair are processed when both overlap the SNP.

**Pitfalls**
- Genotypes do not have to be phased and the sample is always the first one; there is no `-sample` or `-partialVcf` option.
- A read that covers several SNPs writes its calls several times.
- The output files are opened before the VCF index is checked.
- Reads with indels: see "Coordinates of calls".

**Minimal example** (shape from the usage string and `test/unit/allele/SepAlleleMethySCTest.java`; the README has no example)

```
java -Xmx20G -cp "/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.SepAlleleMethySC cell1 phased.vcf.gz cell1.calmd.bam
```

---

## `SepAlleleMethySC2`

**Purpose** — It walks through a read-name-sorted BAM pair by pair; when one end of a pair overlaps a heterozygous SNP the whole pair is assigned to allele A or B and the CpG calls of both ends are written to that allele's file. Choose it over `SepAlleleMethySC` when the mate that does not cover the SNP should be assigned to the allele as well.

**Command**

```
java -Xmx20G -cp "<fat jar>:<ccInference-0.19-slim.jar>:<gatk-package-distribution-3.3.jar>" \
  main.java.edu.mit.compbio.bisulfitehic.utils.SepAlleleMethySC2 [opts] outputPrefix snp_vcf.gz input.sortByReadName.bam
```

Positional arguments (at least 3 tokens are required):
1. `outputPrefix` — prefix of the two output files.
2. `snp_vcf.gz` — SNP VCF.
3. `input.sortByReadName.bam` — alignments sorted by read name.

**Inputs**
- VCF: same rules as `SepAlleleMethySC` (tabix index required, first sample only, heterozygous and not filtered, phasing not required by the code).
- BAM: read sequentially, so the two ends of a pair must be adjacent (read-name sorted); no index is needed. Secondary and supplementary lines are skipped while pairing. Reads need the `MD` tag.
- Chromosome names: the contig name of the read is looked up among the VCF contig names, so both files must use the same names.
- Methylation pattern: fixed reference `CG`, bisulfite filter pattern `CH`; no pattern options.

**Options**

| option | default | what the code does |
|---|---|---|
| `-minBaseQ` | 5 | Minimum base quality of the base at the SNP and of every base used for a call. |
| `-useBadMate` | false | Accepts reads without the proper-pair flag. |
| `-minMapQ` | 30 | Minimum mapping quality, required of both ends. |
| `-turnOffBisulfiteFilter` | false | Skips the bisulfite conversion filter. |
| `-h` | false | No effect (see "Shared behaviour"). |

**Output** — `<outputPrefix>.alleleA_methy.bedgraph.gz` and `<outputPrefix>.alleleB_methy.bedgraph.gz`, gzip-compressed, tab-separated, no header, one line per call in BAM order. Columns are the same as for `SepAlleleMethySC` (contig, position - 1, position of the forward-strand C of the CpG, 0/1); the contig is the contig of the pair.

**Behaviour notes**
- A record is paired with the next primary record; if the read names differ, the first record is dropped and the second starts the next pair.
- Both ends must pass: mapped, primary, `-minMapQ` (30), not vendor-failed, not duplicate, paired, proper pair (unless `-useBadMate`). Both ends must be on the same contig. The distance between the ends is not examined.
- Only one SNP is used per pair: the left-most SNP within the alignment span of the first record, or, when the first record covers none, the left-most SNP within the span of the second record.
- Allele assignment uses the base of the record that covers the SNP, with the same strand rules as `SepAlleleMethySC` (C/T SNP: C->T strand reads skipped; G/A SNP: G->A strand reads skipped; allele A tested first). When the base is missing or below `-minBaseQ` (5) the pair is skipped.
- The calls of both ends go to the allele's file. The call at the SNP position is left out only in the record that covers the SNP.

**Pitfalls**
- The sort order is not checked. In a BAM that is not read-name sorted, every record whose mate is not the next primary record is dropped by the read-name comparison, without any message.
- When the chosen SNP falls in a deletion of the read or has a low-quality base, the pair is skipped; a second SNP on the pair is not tried.
- Pairs with the ends on different contigs are always skipped.
- Overlapping ends write a call for the same CpG twice.
- Genotypes do not have to be phased and the sample is always the first one.
- Reads with indels: see "Coordinates of calls".

**Minimal example** (shape from the usage string and `test/unit/allele/SepAlleleMethySC2Test.java`; the README has no example)

```
java -Xmx20G -cp "/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.SepAlleleMethySC2 cell1 phased.vcf.gz cell1.sortByReadName.calmd.bam
```

---

## `EnhanceHicMatrix`

**Purpose** — For pairs of genomic windows on the same chromosome it takes the read pairs that have one end in each window, computes the methylation level of each end, and reports the Pearson correlation between the two ends across read pairs, together with a permuted control. It does not use a VCF and is not allele-specific; choose it to score methylation concordance between two anchors (from a bedpe file) or between all bins of a chromosome.

**Command**

```
java -Xmx20G -cp "<fat jar>:<ccInference-0.19-slim.jar>:<gatk-package-distribution-3.3.jar>" \
  main.java.edu.mit.compbio.bisulfitehic.utils.EnhanceHicMatrix [opts] outputPrefix hg19.fa.fai cpg_list.bed[.gz] input.bam
```

Positional arguments (at least 4 tokens are required):
1. `outputPrefix` — prefix of the output file.
2. `hg19.fa.fai` — chromosome sizes.
3. `cpg_list.bed[.gz]` — intervals used to select reads.
4. `input.bam` — alignments.

**Inputs**
- Chromosome size file: tab-separated, column 1 = contig, column 2 = length (a `.fai` works); `#` lines and lines with fewer than 2 columns are skipped; gzip when the name ends in `.gz`. It is always read but only used to generate windows when `-bedpe` is not given.
- CpG list: tab-separated, at least 3 columns (contig, start, end); `#` lines and shorter lines are skipped; gzip when the name ends in `.gz`. It works as a gate: a read is used only when its alignment span overlaps at least one listed interval on its contig. Start and end are used as written (no 0-based/1-based conversion). Column 6 (strand) is parsed, but the stored value is never read. The methylation calls do not come from this file; they come from the pattern search on the read.
- BAM: coordinate-sorted and indexed (region queries); reads need the `MD` tag; the two ends of a pair must have the same read name.
- `-bedpe` file: tab-separated, columns 1-6 = contig1, start1, end1, contig2, start2, end2; `#` lines are skipped; gzip when the name ends in `.gz`.
- Chromosome names: the same names must be used in the bedpe/size file, the CpG list and the BAM.
- Pattern options: see "Shared behaviour".

**Options**

| option | default | what the code does |
|---|---|---|
| `-minBaseQ` | 5 | Minimum base quality of a base used for a call. |
| `-turnOffBisulfiteFilter` | false | Skips the bisulfite conversion filter. |
| `-useBadMate` | false | Accepts reads without the proper-pair flag. |
| `-useUnpaired` | false | Accepts reads without the paired flag (treated as first-of-pair). See Pitfalls: it needs `-useBadMate`. |
| `-minMapQ` | 30 | Minimum mapping quality. |
| `-minDistEnds` | 20000 | A window pair is dropped when `abs(start1-end2)` or `abs(end1-start2)` is below this value. It is applied to window/anchor coordinates, not to reads. |
| `-maxDistEnds` | 2000000 | A window pair is dropped when `abs(start1-start2)` or `abs(end1-end2)` is above this value. |
| `-resolution` | 25000 | Bin size of the generated window pairs. Used only without `-bedpe`. |
| `-window` | 5000 | Size of the sub-windows that scan each anchor. Used only with `-bedpe`. |
| `-slidingWindow` | false | Sub-window step is `window/10` (integer division) instead of `window`. Used only with `-bedpe`. |
| `-minReadsInWindow` | 3 | Minimum number of reads in each window, and minimum number of read pairs shared by the two windows. |
| `-bedpe` | null | File of anchor pairs to analyse. Without it all bin pairs of every chromosome in the size file are generated. |
| `-chrom` | null | Keeps only this contig (case-insensitive) in the CpG list, the bedpe rows and the generated windows. |
| `-minCgNum` | 1 | Minimum number of calls in a read for the read to be used. |
| `-permutations` | 10 | Number of random re-pairings concatenated for the permuted correlation. |
| `-skipProcessFirstRow` | false | Skips the first line of the bedpe file. |
| `-methyPatternSearch` | `CG` | Reference pattern for methylation calls. |
| `-methyPatternCytPos` | 1 | 1-based position of the inspected base in the pattern. |
| `-bsConvPatternSearch` | `CH` | Pattern used by the bisulfite conversion filter. |
| `-bsConvPatternCytPos` | 1 | Passed to the filter, but it does not change the filter result (see "Shared behaviour"). |
| `-h` | false | No effect (see "Shared behaviour"). |

**Output** — one file, `<outputPrefix>.enhancedHicMatrix.txt.gz`, gzip-compressed, tab-separated, no header. A row is written only when the observed correlation (column 8) is a number. Numbers are printed as Java doubles (for example `4.0`, `NaN`).

| col | meaning |
|---|---|
| 1 | contig of window 1 |
| 2 | start of window 1 |
| 3 | end of window 1 |
| 4 | contig of window 2 |
| 5 | start of window 2 |
| 6 | end of window 2 |
| 7 | name: `contig:start1:end1:start2:end2`. Without `-bedpe` these are the window coordinates. With `-bedpe` they are the coordinates of the bedpe anchors, while columns 2, 3, 5, 6 are the sub-windows. |
| 8 | Pearson correlation between the methylation level of the end in window 1 and of the end in window 2, over the read names present in both windows |
| 9 | p-value of column 8 |
| 10 | number of read names present in both windows |
| 11 | Pearson correlation after random re-pairing: in each of `-permutations` rounds every read of window 1 is paired with a randomly drawn read of window 2; all rounds are pooled into one correlation |
| 12 | p-value of column 11 |
| 13 | number of points in the pooled permutation (`permutations` x number of reads in window 1) |

Columns 11-13 are `NaN`, `NaN`, `0.0` when the p-value of column 8 is `NaN`.

**Behaviour notes**
- For each window pair the BAM is queried twice (window 1, window 2). A read is used when its span overlaps a CpG-list interval and the helper accepts it: mapped, primary, not duplicate, not vendor-failed, paired (unless `-useUnpaired`), proper pair (unless `-useBadMate`), mapping quality >= `-minMapQ` (30), and not rejected by the bisulfite filter.
- The methylation level of a read is methylated calls / all calls over the whole read (not only the part inside the window, and not only the listed CpGs). A read with fewer than `-minCgNum` calls is dropped.
- The distance between the two ends of a pair is never measured; a pair contributes when its two ends fall in the two windows.
- Without `-bedpe`: for every chromosome, bins of `-resolution` starting at 0 are paired with every later bin, subject to `-minDistEnds`/`-maxDistEnds`; each bin pair is tested once as a whole.
- With `-bedpe`: rows are kept when both anchors are on the same contig and pass `-minDistEnds`/`-maxDistEnds`. Each anchor is scanned with sub-windows of `-window`, and every sub-window pair is tested.
- The random generator is not seeded, so columns 11-13 differ between runs.

**Pitfalls**
- `-window` and `-slidingWindow` have no effect without `-bedpe`, and `-resolution` has no effect with `-bedpe`.
- With `-bedpe`, an anchor shorter than `-window` produces no sub-window and therefore no output for that row.
- With `-slidingWindow` and `-window` below 10 the step is 0 and the scan never advances.
- `-useUnpaired` without `-useBadMate` makes htsjdk throw `IllegalStateException: Inappropriate call if not paired read` on an unpaired read (the helper asks for the proper-pair flag; htsjdk `SAMRecord.getProperPairFlag`).
- `-skipProcessFirstRow` skips line 1 only; if line 1 is a `#` comment, a header on line 2 is not skipped.
- Without `-bedpe` all window pairs of the selected chromosomes are built in memory before the BAM is read.
- The class has a private `failFlagFilter` method that is never called; read filtering happens in the helper.

**Minimal example** (shape from the usage string and `test/unit/allele/EnhanceHicMatrixTest.java`; the README has no example)

```
java -Xmx20G -cp "/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.EnhanceHicMatrix sample genome.fa.fai cpg_list.bed.gz sample.calmd.bam -bedpe loops.bedpe -useBadMate
```

---

## Known limitations

Found while reading the source for this guide (November 2026). They are described here, not fixed.

All five tools

1. `-h` does nothing. Usage is printed only when too few tokens are given; the count includes options, so enough tokens with too few positional arguments ends in an index error. Usage errors exit with status 0.
2. Insertions and deletions shift the coordinates of methylation calls after them: coordinates are alignment start plus offset in the clipped read.
3. `-useUnpaired` cannot add reads (`LongRangeAsm`), or needs `-useBadMate` (`EnhanceHicMatrix`); otherwise an unpaired read stops the run with `IllegalStateException`.
4. Not verified: behaviour when the VCF header has no `##contig` lines (the sequence dictionary is taken from the VCF), and when a contig is missing from the VCF header or the BAM header.

`LongRangeAsm`

5. "ref" and "alt" mean allele A and allele B, the first and second allele of the phased genotype, not REF and ALT. This applies to `-coverageRef`, `-coverageAlt` and to the `output.ref` / `output.alt` files. The coverage thresholds are applied to the `*_local_count` columns.
6. `-qual` is parsed and never used (also in `AlleImbalanceExpr`). `-bsConvPatternCytPos` does not change the result (also in `EnhanceHicMatrix`).
7. The usage text of `-methyPattern2ndSearch` says "default: CG"; the real default is off. The usage text of `-minBaseQ` says "mapping quality"; the default is 20 here and 5 in the other tools.
8. The per-cytosine tables (`output.ref`, `output.alt`) are filled before the coverage test and once per overlapping SNP. They include SNPs absent from the main table and can count a read several times.
9. `End2_anchor_point` and `MaxEnd2Dist` include reads that match neither allele. With `-diffChr` they mix coordinates from different contigs, and the anchor list carries no contig name.
10. A methylation call at the SNP position itself is kept (the two `SepAlleleMethySC` tools drop it).
11. On G-to-A strand reads the reported coordinate is the G of the CpG on the forward strand, so the two strands of one CpG land on adjacent positions. The two `SepAlleleMethySC` tools report the C for both strands.

`SepAlleleMethySC`, `SepAlleleMethySC2`

12. Phasing is not required: the code tests only that the genotype is heterozygous. The sample is fixed to the first one in the VCF.
13. On G-to-A strand reads a SNP at the G of a CpG is not skipped, because the skip compares against the reported C offset.

`AlleImbalanceExpr`

14. Duplicate-flagged reads are not dropped (the other allele tools drop them).
15. `-consistentPrefix` only adds `chr` to the contig name used in the BAM query.

`EnhanceHicMatrix`

16. `-window` and `-slidingWindow` work only with `-bedpe`; `-resolution` works only without it. `-slidingWindow` with `-window` below 10 gives a step of 0 and the run does not end.
17. The strand column of the CpG list is read and not used. Read methylation is taken over the whole read, not over the window or the listed CpGs.
18. `-skipProcessFirstRow` skips only the first physical line of the BEDPE file.

Not covered by an example or test

19. The README has an example only for `LongRangeAsm`. The examples for the other four tools on this page are built from the usage strings and the unit tests.
