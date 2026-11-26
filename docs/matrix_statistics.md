# Matrix and statistics tools

Part of the [usage guide](USAGE.md). These Java tools build CpG-by-CpG or window-by-window concordance tables and test them.

| Tool | Input | Output | Choose it when |
|---|---|---|---|
| [`SingleBaseResHicMatrix`](#singlebasereshicmatrix) | CpG list + BAM | one row per CpG pair linked by read pairs, with concordance statistics (17 columns) | default choice; loads the reads of a chromosome once |
| [`SingleBaseResHicMatrixV2`](#singlebasereshicmatrixv2) | CpG list + BAM | same 17 columns, same file name | memory is short; queries the BAM again for every CpG |
| [`SingleCellTest`](#singlecelltest) | cell-by-site table | site-by-site score matrix | per-cell methylation table, correlation or mutual information |
| [`TestIndepPoiBinom`](#testindeppoibinom) | cell-by-site table | site-by-site matrix | observed against expected co-methylation, optional permutation |
| [`TestIndepPoiBinomOnBam`](#testindeppoibinomonbam) | SNP VCF + per-cell name-sorted BAMs | window-by-window matrix | the same test from BAMs, split by allele |

- `PoissonBinomial2` and `FisherExactTest` are library classes, not tools.
- All five tools are covered by `test/unit/stats/`.
- These tools are fixed to CpG (reference `C` followed by `G`) with a `CH` conversion filter; they have no pattern options.
- Sentences that start with "Not verified:" rest on library behaviour that was not checked.

---

## Common to all five tools

### Classpath

The main classes are in the package `main.java.edu.mit.compbio.bisulfitehic.utils`. The README runs the `utils` tools with this classpath shape:

```
java -Xmx20G -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/ccInference-0.19-slim.jar:/path/to/bisulfiteHic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.<ClassName> [options] <positional arguments>
```

The jar `bisulfitehic-0.38-jar-with-dependencies.jar` contains htsjdk, args4j, commons-math3 and log4j, but it does not contain the classes of `lib/ccInference-0.19-slim.jar`, `lib/gatk-package-distribution-3.3.jar`, `lib/JavaMI.jar` or `lib/java-genomics-io.jar` (these are `system`-scope dependencies in `pom.xml`; a listing of the jar shows no such classes). Which extra jar each tool refers to in its source:

| Tool | Extra jars referenced by the source |
|---|---|
| `SingleBaseResHicMatrix`, `SingleBaseResHicMatrixV2` | `ccInference-0.19-slim.jar` and `gatk-package-distribution-3.3.jar` (through `BisulfiteHicUtils`) |
| `SingleCellTest` | `JavaMI.jar` (`JavaMI.MutualInformation`, called for `-scoreMode 2`) |
| `TestIndepPoiBinom` | none directly; with `-permutation` > 0 it calls `BisulfiteHicUtils.findPosition`, and `BisulfiteHicUtils` refers to ccInference and GATK classes |
| `TestIndepPoiBinomOnBam` | `ccInference-0.19-slim.jar` and `gatk-package-distribution-3.3.jar` |

The unit-test runner (`test/unit/stats/run.sh`) uses one classpath for all of them: the fat jar first, then `lib/*`. The same shape works as a single recipe:

```
java -Xmx20G -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/*" main.java.edu.mit.compbio.bisulfitehic.utils.<ClassName> ...
```

Not verified: whether `SingleCellTest` with `-scoreMode 0` or `1`, or `TestIndepPoiBinom` without `-permutation`, start without the extra jars on the classpath (this depends on JVM class resolution and was not tested).

Not verified: whether the order of the jars matters (the README and `run.sh` both put the fat jar first; the GATK jar may bundle other versions of the same libraries).

The README lists Oracle JDK 8 as the prerequisite. The unit tests compile with `--release 8` and run on OpenJDK 11.

### Help and argument checking

In all five classes the test `if(help || args.length < N)` runs before the command line is parsed, so the `-h` field is still `false` at that point. Usage (the `USAGE` string and the option list) is printed to standard error only when fewer than `N` command-line tokens are given (`N` = 3 for the two `SingleBaseResHicMatrix` classes, 2 for the other three), or when args4j rejects an option. In these cases the method returns normally and no output file is created. `args.length` counts option tokens too, so a command line with enough option tokens but missing positional arguments passes this check and then fails at `arguments.get(i)`. Giving `-h` together with enough tokens does not print help; the tool runs. Extra positional arguments are ignored.

### Output number format

All numeric fields are written with Java's default `double` to string conversion. Counts therefore appear as `3.0`, undefined values as `NaN`, a division of a positive number by zero as `Infinity`, and very small or large values in scientific notation such as `1.0E-4`.

### Read filtering and CpG methylation calls (`BisulfiteHicUtils`)

`SingleBaseResHicMatrix`, `SingleBaseResHicMatrixV2` and `TestIndepPoiBinomOnBam` obtain per-read CpG calls from `BisulfiteHicUtils.readsMethySummaryAtEachLoc(read, minMapQ, minBaseQ, turnOffBisulfiteFilter, useBadMate)`. As implemented:

- It returns `null` for a read that is unmapped, secondary (flag 0x100), duplicate, QC-failed, not paired, has mapping quality below `minMapQ`, or lacks the proper-pair flag while `useBadMate` is false. The supplementary flag (0x800) is not tested.
- The reference bases under the read are rebuilt from the read itself with htsjdk `SequenceUtil.makeReferenceFromAlignment(read, false)`. No reference FASTA is opened by these tools. The test data carry the reference only in the `MD` tag. Not verified: a read without an `MD` tag makes htsjdk throw (the bundled htsjdk class contains the message "Cannot create reference from SAMRecord with no MD tag"). The README adds `MD` with `samtools calmd`.
- Soft-clipped bases are removed; deleted reference bases are skipped; inserted read bases stay in the arrays. If the rebuilt reference and the clipped read differ in length an `IllegalArgumentException` is thrown.
- Unless `turnOffBisulfiteFilter` is set, `CcInferenceUtils.bisulfiteIncompleteReads(readBases, refBases, negStrand, secondPair, false, true, true, 0.1, 0.4, 0, "CH", 1)` is called. A negative return value makes the read contribute no CpG at all; a non-negative value is used as the offset from which CpG calls are reported (for G-to-A reads the last `offset` bases are excluded instead). Not verified: what that function computes; it lives in `lib/ccInference-0.19-slim.jar`, whose source is not in this repository.
- The pattern is fixed: a reference `C` followed by a reference `G`. These five tools have no `-methyPatternSearch` / `-bsConvPatternSearch` options.
- Reads that are first-of-pair on the forward strand or second-of-pair on the reverse strand are read as C-to-T reads: read base `C` at the reference `C` gives value 1 (methylated), read base `T` gives 0. All other reads are read as G-to-A reads: read base `G` at the reference `G` of the CpG gives 1, read base `A` gives 0. Any other base, or a base quality below `minBaseQ`, gives no call.
- The result is a map from the 0-based offset of the CpG's `C` (counted in aligned read bases from the first aligned base, forward reference strand, for both read types) to the value 0 or 1.

`BisulfiteHicUtils.findPosition(sortedList, observed, 0, n-1)` is used for all permutation ranks. On an ascending list of `n` values it returns 0 when `observed` is less than or equal to the smallest value, `n-1` when `observed` is greater than the largest value, and otherwise an index `k` with `list[k-1] <= observed <= list[k]` found by bisection. The returned rank is therefore never `n`.

---

## `SingleBaseResHicMatrix`

**Purpose** — For every pair of CpG positions on the same chromosome it counts the read names that have a methylation call at both positions, and reports how often the two calls agree, their correlation and entropy, and the same quantities under random re-pairing of reads. It works on one indexed BAM and has no notion of cells.

**Command**

```
java -Xmx20G -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/ccInference-0.19-slim.jar:/path/to/bisulfiteHic/lib/gatk-package-distribution-3.3.jar" \
  main.java.edu.mit.compbio.bisulfitehic.utils.SingleBaseResHicMatrix [options] outputPrefix cpg_list.bed[.gz] input.bam
```

Positional arguments, in order:

1. `outputPrefix` — prefix of the output file (`arguments.get(0)`). The directory must already exist.
2. `cpg_list.bed[.gz]` — CpG position list (`arguments.get(1)`). Read through gzip when the name ends with `.gz`.
3. `input.bam` — alignment file (`arguments.get(2)`).

**Inputs**

- `input.bam`: opened with htsjdk with validation stringency SILENT and queried by region (`queryOverlapping`), one multi-interval query per chromosome. The unit test uses a coordinate-sorted BAM with a `.bai` index. Not verified: a BAM without an index fails (the bundled htsjdk contains the message "No index is available for this BAM file."). Reads need the paired flag, the proper-pair flag and (see the common section) an `MD` tag. No other tag is read. Read groups are ignored. The two ends of a contact are linked by identical read name.
- `cpg_list.bed[.gz]`: tab-separated. Lines starting with `#` and lines with fewer than 3 columns are skipped. There is no first-row skip for this file. Columns used (0-based index in the code): 0 = chromosome, 1 = start, 2 = end. Column index 5 (strand `+`/`-`) is parsed and stored but never used afterwards.
  - Coordinate convention: the start column is used directly as the 1-based reference position at which the methylation call is looked up, and the calls are keyed by the position of the CpG's `C`. The start must therefore be the 1-based position of the `C`. The unit test writes each CpG as `chr1<TAB>C_position<TAB>C_position+1`. With a standard 0-based BED start (one less than the 1-based `C` position) the lookup hits the base before the `C`, where no call exists, so no CpG is matched. The test README states that support for a standard 0-based BED is an open question.
  - The BAM query interval for a CpG is `start+1 .. end`, so `end` must be at least `start+1`.
  - Two lines with identical start and end collapse into one CpG.
  - The chromosome name must be a sequence name in the BAM header (`dict.getSequenceIndex(chr)`). Not verified: what happens for a name missing from the BAM header (the index is then -1 and is passed to the query).
- `-bedpe` file (optional, `.gz` when the name ends with `.gz`): tab-separated; lines starting with `#` are skipped; with `-skipProcessFirstRow` the first non-`#` line is skipped. Columns (0-based): 0 chr1, 1 start1, 2 end1, 3 chr2, 4 start2, 5 end2. Coordinates are used as written. A row is dropped when chr1 and chr2 differ (compared ignoring case), or when `|start1-start2|` or `|end1-end2|` is below `-minDistEnds` or above `-maxDistEnds`. Without `-skipProcessFirstRow`, a header line that does not start with `#` causes a `NumberFormatException`.
- `-regions` file (optional, `.gz` when the name ends with `.gz`): tab-separated BED; lines starting with `#` are skipped; with `-skipProcessFirstRow` the first non-`#` line is skipped as well. Columns (0-based): 0 chr, 1 start, 2 end. Coordinates are used as written.
- Chromosome names in `-bedpe` and `-regions` must be spelled exactly like those in the CpG list (hash-map lookup, case-sensitive).
- Methylation pattern: fixed CpG; see the common section. There are no pattern options.

**Options**

| Option | Default | What it does in the code |
|---|---|---|
| `-minBaseQ` | 5 | Minimum base quality of the read base at the CpG for a call to be made. |
| `-turnOffBisulfiteFilter` | false | Skips the call to `CcInferenceUtils.bisulfiteIncompleteReads`; without it a read for which that function returns a negative value contributes no call. |
| `-useBadMate` | false | Parsed and passed to `readsMethySummaryAtEachLoc`, but it has no effect: the class's own `failFlagFilter` discards every read without the proper-pair flag before that call. |
| `-minMapQ` | 30 | Reads with mapping quality below this are discarded. |
| `-minDistEnds` | 1000 | A CpG pair is skipped when `|end_i - end_j|` is below this. Also filters `-bedpe` rows (see Inputs). |
| `-maxDistEnds` | 200000000 | A CpG pair is skipped when `|end_i - end_j|` is above this. Also filters `-bedpe` rows, and with `-regions` it is the margin around each CpG within which a region must lie for the CpG to be loaded. |
| `-bedpe` | null | Restricts the output to CpG pairs where CpG i overlaps one anchor of a bedpe row and CpG j intersects a partner anchor of an anchor overlapping CpG i. Only chromosomes present in the kept bedpe rows are processed. The `name` column then holds the anchor coordinates. |
| `-minReadsInWindow` | 3 | (a) A CpG is kept only if at least this many read names have a call at it. (b) A CpG pair is scored only if at least this many read names (and at least one) have a call at both CpGs. |
| `-bgBam` | null | Parsed; the file is opened and closed, so it must be readable. It has no effect on the output: the code that would fill the background table loops over a list that was cleared just before, so the table stays empty and the background branch is never taken. |
| `-permutations` | 100 | Number of resampling rounds per scored CpG pair. With 0, `permutated_ratio`, `cor_random`, `pvalue_random`, `entropy_random` and `entropy_permuation_p_value` are `NaN`. |
| `-skipProcessFirstRow` | false | Skips the first non-`#` line of the `-bedpe` file and of the `-regions` file. It does not apply to the CpG list. |
| `-regions` | null | (a) A CpG is loaded only if a region on its chromosome overlaps `[start - maxDistEnds, end + maxDistEnds]`. (b) A CpG pair is kept only if at least one of the two CpGs overlaps a region. |
| `-h` | false | See "Help and argument checking": the flag itself is never effective; usage is printed when fewer than 3 tokens are given. |

**Output**

One file: `<outputPrefix>.singleBaseResHicMatrix.bed.gz`, gzip-compressed, tab-separated, UTF-8. The first line is a header that starts with `#`. One row is written for each CpG pair (i, j) on the same chromosome, with CpG i before CpG j in position order, that passes the filters and has at least `max(1, -minReadsInWindow)` shared read names. Pairs with fewer shared read names are not written. Chromosomes appear in hash-map order, not sorted.

| # | Header label | Meaning | Missing value |
|---|---|---|---|
| 1 | `#chr` | Chromosome. | never |
| 2 | `start1` | Start of CpG i as given in the CpG list. | never |
| 3 | `end1` | End of CpG i as given in the CpG list. | never |
| 4 | `start2` | Start of CpG j. | never |
| 5 | `end2` | End of CpG j. | never |
| 6 | `name` | `chr:start1:end1:start2:end2`. With `-bedpe`: `chr:a1start:a1end:a2start:a2end`, where a1 is the first bedpe anchor overlapping CpG i and a2 is the first partner anchor intersecting CpG j. | never |
| 7 | `concordCount` | Number of shared read names whose calls at the two CpGs are equal. Written as a double (`3.0`). | never |
| 8 | `disconcordCount` | Number of shared read names whose calls differ. Written as a double. | never |
| 9 | `observed_ratio` | `log10(concordCount + 0.0001) - log10(disconcordCount + 0.0001)`. | never in a written row |
| 10 | `permutated_ratio` | Mean of the same ratio over the permutations. | `NaN` |
| 11 | `correlation` | Pearson correlation of the 0/1 calls at CpG i and CpG j over the shared read names (commons-math `PearsonsCorrelation`). | `NaN` (for example when all calls at one CpG are equal) |
| 12 | `correlation_pvalue` | Entry (0,1) of `PearsonsCorrelation.getCorrelationPValues()` for that correlation; only computed when the correlation is not NaN. | `NaN` |
| 13 | `cor_random` | Mean of the non-NaN Pearson correlations of the permutations. | `NaN` |
| 14 | `pvalue_random` | Mean of the correlation p-values of the permutations. | `NaN` |
| 15 | `entropy` | Shannon entropy (base 2) of the joint call patterns (`00`, `01`, `10`, `11`) over the shared read names, divided by `log2(number of shared read names)`. | `NaN` when exactly one read name is shared (division by `log2(1)`), which requires `-minReadsInWindow 1` |
| 16 | `entropy_random` | Mean of the same normalised entropy over the permutations. | `NaN` |
| 17 | `entropy_permuation_p_value` | `findPosition` rank of the observed entropy in the ascending list of permutation entropies, divided by the number of permutations. | `NaN` when `-permutations 0` |

The header spells the labels exactly as shown (`disconcordCount`, `permutated_ratio`, `entropy_permuation_p_value`).

**Method**

For each chromosome, all reads overlapping any listed CpG are fetched once and, per CpG, a table `read name -> call (0/1)` is built. For a CpG pair, the shared read names are those present in both tables. Concordant means equal calls. The Pearson correlation and its p-value come from commons-math. The entropy is the normalised joint-pattern entropy described above.

Null model: for each of `-permutations` rounds, `T` draws are made, where `T` is the number of shared read names. Each draw takes one read name uniformly at random (with replacement) from all read names with a call at CpG i, and independently one from all read names with a call at CpG j, and pairs their calls. The ratio, correlation, correlation p-value and entropy are computed on these `T` random pairs. The random number generator is `new MersenneTwister()` without a seed; there is no seed option, so the permutation columns differ between runs.

Not verified: the exact definition of the commons-math correlation p-value (expected: two-sided test based on the t distribution).

**Behaviour notes**

- A read is used only if it is mapped, primary (flag 0x100 not set), not a duplicate, not QC-failed, paired, has the proper-pair flag and mapping quality at least `-minMapQ`. Supplementary alignments are not excluded.
- If the same read name gives a call at the same CpG more than once (for example overlapping mates), the first record seen is kept.
- Calls from both read types (C-to-T and G-to-A) are stored at the `C` position, so the two strands of a CpG are combined.
- Only intra-chromosomal CpG pairs are produced. A shared read name can come from the two mates of a pair or from one read covering both CpGs.
- Memory: for one chromosome at a time, the tool holds a read-name table for every CpG that has calls. Run time: every pair of kept CpGs on a chromosome within the distance limits is visited (the default limits, 1 kb to 200 Mb, allow nearly all pairs); permutations are run only for pairs with enough shared read names.
- Progress is logged every 1,000,000 visited CpG pairs, and the writer is flushed at that point.
- The output file and its header are written even if no pair qualifies.

**Pitfalls**

- The CpG list must use the 1-based position of the `C` as the start column (see Inputs). A standard 0-based BED produces a file that contains only the header.
- `-useBadMate` and `-bgBam` do not change the result.
- Reads without the proper-pair flag are always dropped. Not verified: whether the long-range pairs in a given BAM carry that flag depends on how the BAM was produced; this was not checked here.
- `-skipProcessFirstRow` acts on both the `-bedpe` and the `-regions` file; a header-less `-regions` file loses its first region when the option is given for the sake of the bedpe header.
- With `-bedpe`, CpGs on a chromosome named in the bedpe are still all loaded (the code tests `overlappers(...) == null`, which the code's own comments say never happens); the restriction to anchors is applied later in the pair loop, so `-bedpe` reduces output but not the read loading for those chromosomes.
- The bedpe restriction is applied only when at least one bedpe row survives the row filter. If the file is empty, or every row is dropped (different chromosomes, or anchor distance outside `-minDistEnds`/`-maxDistEnds`), the run proceeds without any restriction and without a warning. The same holds for a `-regions` file that yields no region.
- `SingleBaseResHicMatrixV2` writes to the same file name; using the same prefix for both overwrites the file.
- Not verified: the multi-interval BAM query may be rejected by htsjdk when two CpG query intervals (`start+1 .. end`) overlap or abut (the bundled htsjdk contains the messages "List of intervals is not optimized: %s intersects %s" and "... abuts %s"). The test data have CpGs 4 kb apart and do not exercise this.
- The permutation columns are not reproducible between runs.

**Minimal example** (invocation of `test/unit/stats/SingleBaseResHicMatrixTest.java`, case 1; `-turnOffBisulfiteFilter` is part of the test invocation and is optional)

```
java -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/*" \
  main.java.edu.mit.compbio.bisulfitehic.utils.SingleBaseResHicMatrix -turnOffBisulfiteFilter /path/to/out/prefix /path/to/cpg.bed /path/to/reads.bam
```

This writes `/path/to/out/prefix.singleBaseResHicMatrix.bed.gz`. In the test, the CpG list is `chr1<TAB>1000<TAB>1001` and so on, and the CpG pair 1000/5000 with three concordant read pairs gives a row whose columns 7 and 8 are `3.0` and `0.0`.

---

## `SingleBaseResHicMatrixV2`

**Purpose** — It answers the same question as `SingleBaseResHicMatrix` and writes the same 17 columns to the same file name. The difference is how reads are obtained: V2 keeps nothing in memory and queries the BAM again for every CpG of every pair, whereas V1 loads all reads of a chromosome once. V2 has no `-regions` option and does not pre-filter single CpGs by read count. Based on these code differences, V1 is the one to use unless memory is the limiting factor and the CpG list is small or restricted with `-bedpe`; V2 performs one BAM query per CpG pair.

**Command**

```
java -Xmx20G -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/ccInference-0.19-slim.jar:/path/to/bisulfiteHic/lib/gatk-package-distribution-3.3.jar" \
  main.java.edu.mit.compbio.bisulfitehic.utils.SingleBaseResHicMatrixV2 [options] outputPrefix cpg_list.bed[.gz] input.bam
```

Positional arguments, in order: `outputPrefix` (`arguments.get(0)`), `cpg_list.bed[.gz]` (`arguments.get(1)`), `input.bam` (`arguments.get(2)`); same meaning as in V1. The class name ends in a capital `V2`; the usage string prints `SingleBaseResHicMatrixv2`.

**Inputs**

- `input.bam`: as for V1 (region queries, so an indexed coordinate-sorted BAM as in the unit test; paired and proper-pair flags and `MD` tag needed; read groups ignored; ends linked by read name). V2 queries by chromosome name: `queryOverlapping(chr, start, end)` with the CpG list's start and end as written.
- `cpg_list.bed[.gz]`: same parsing as V1 (skip `#` lines and lines with fewer than 3 columns; 0-based column indexes 0, 1, 2; strand column index 5 parsed but unused; no first-row skip). The same coordinate convention applies: the start column must be the 1-based position of the CpG's `C`.
- `-bedpe` file: same format, same row filter and same `-skipProcessFirstRow` handling as V1.
- There is no `-regions` input.
- Methylation pattern: fixed CpG; no pattern options.

**Options**

| Option | Default | What it does in the code |
|---|---|---|
| `-minBaseQ` | 5 | Minimum base quality of the read base at the CpG for a call to be made. |
| `-turnOffBisulfiteFilter` | false | Skips the call to `CcInferenceUtils.bisulfiteIncompleteReads`. |
| `-useBadMate` | false | Parsed and passed on, but without effect: `failFlagFilter` discards every read without the proper-pair flag first. |
| `-minMapQ` | 30 | Reads with mapping quality below this are discarded. |
| `-minDistEnds` | 1000 | A CpG pair is skipped when `|end_i - end_j|` is below this. Also filters `-bedpe` rows. |
| `-maxDistEnds` | 2000000000 | A CpG pair is skipped when `|end_i - end_j|` is above this. Also filters `-bedpe` rows. (The V1 default is 200000000.) |
| `-bedpe` | null | Same restriction and same effect on the `name` column as in V1. |
| `-bgBam` | null | Parsed only. The field is never read anywhere in the class; the file is not even opened. |
| `-minReadsInWindow` | 3 | A CpG pair is scored only if at least this many read names (and at least one) have a call at both CpGs. There is no per-CpG pre-filter. |
| `-permutations` | 100 | Number of resampling rounds per scored pair; with 0 the permutation columns are `NaN`. |
| `-skipProcessFirstRow` | false | Skips the first non-`#` line of the `-bedpe` file. It does not apply to the CpG list. |
| `-h` | false | See "Help and argument checking": the flag itself is never effective; usage is printed when fewer than 3 tokens are given. |

**Output**

One file: `<outputPrefix>.singleBaseResHicMatrix.bed.gz` (the same suffix as V1), gzip-compressed, tab-separated, with the same `#chr ...` header line and the same 17 columns as described for `SingleBaseResHicMatrix`. A row is written only when `observed_ratio` is not NaN, that is when at least `max(1, -minReadsInWindow)` read names are shared. Chromosomes appear in hash-map order.

| # | Header label | Meaning | Missing value |
|---|---|---|---|
| 1 | `#chr` | Chromosome. | never |
| 2 | `start1` | Start of CpG i as given in the CpG list. | never |
| 3 | `end1` | End of CpG i. | never |
| 4 | `start2` | Start of CpG j. | never |
| 5 | `end2` | End of CpG j. | never |
| 6 | `name` | `chr:start1:end1:start2:end2`, or with `-bedpe` the coordinates of the two anchors, `chr:a1start:a1end:a2start:a2end`. | never |
| 7 | `concordCount` | Shared read names with equal calls (written as a double). | never |
| 8 | `disconcordCount` | Shared read names with different calls (written as a double). | never |
| 9 | `observed_ratio` | `log10(concordCount + 0.0001) - log10(disconcordCount + 0.0001)`. | never in a written row |
| 10 | `permutated_ratio` | Mean of that ratio over the permutations. | `NaN` |
| 11 | `correlation` | Pearson correlation of the 0/1 calls over the shared read names. | `NaN` |
| 12 | `correlation_pvalue` | p-value of that correlation from commons-math. | `NaN` |
| 13 | `cor_random` | Mean of the non-NaN permutation correlations. | `NaN` |
| 14 | `pvalue_random` | Mean of the permutation correlation p-values. | `NaN` |
| 15 | `entropy` | Normalised entropy of the joint call patterns over the shared read names. | `NaN` when exactly one read name is shared |
| 16 | `entropy_random` | Mean of the permutation entropies. | `NaN` |
| 17 | `entropy_permuation_p_value` | `findPosition` rank of the observed entropy among the sorted permutation entropies, divided by their number. | `NaN` when `-permutations 0` |

**Method**

Identical formulas to V1 (the `cpgContactFreq` method is the same code): shared read names, concordant and discordant counts, log10 ratio with a 0.0001 pseudo-count, Pearson correlation and p-value, normalised joint-pattern entropy, and the same resampling null model with an unseeded `MersenneTwister`.

**Behaviour notes**

- Same read filter as V1 (mapped, primary, not duplicate, not QC-failed, paired, proper pair, mapping quality at least `-minMapQ`); first record per read name wins.
- For each CpG i the BAM is queried once, and for each partner CpG j that passes the distance and bedpe checks the BAM is queried again. Nothing is cached, so the number of BAM queries grows with the number of CpG pairs. Memory use is small.
- A list of query intervals and the sequence-dictionary index of each chromosome are computed but not used.
- The end-of-run log line says "SingleBaseResHicMatrix's running time".

**Pitfalls**

- Same CpG-list coordinate convention as V1; a standard 0-based BED gives only the header.
- `-useBadMate` and `-bgBam` have no effect.
- Same output file name as V1.
- Run time grows with the square of the number of CpGs per chromosome, each step being a BAM query.
- The permutation columns are not reproducible between runs.

**Minimal example** (invocation of `test/unit/stats/SingleBaseResHicMatrixTest.java` with argument `v2`, case 3)

```
java -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/*" \
  main.java.edu.mit.compbio.bisulfitehic.utils.SingleBaseResHicMatrixV2 -turnOffBisulfiteFilter -skipProcessFirstRow -bedpe /path/to/loops.bedpe /path/to/out/prefix /path/to/cpg.bed /path/to/reads.bam
```

In the test, `loops.bedpe` has one header line and the row `chr1 900 1100 chr1 8900 9100`; the only output row is the CpG pair 1000/9000 and its `name` column is `chr1:900:1100:8900:9100`.

---

## `SingleCellTest`

**Purpose** — From a text matrix with one row per CpG (or region) and one value column per cell, it computes a score for every pair of rows across the cells that have a value in both rows, and writes the full square matrix. The score is a Pearson correlation, a variance ratio or a mutual information, chosen with `-scoreMode`.

**Command**

```
java -Xmx20G -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/JavaMI.jar" \
  main.java.edu.mit.compbio.bisulfitehic.utils.SingleCellTest [options] input.bed.gz outputPrefix
```

Positional arguments, in order:

1. `input.bed.gz` — the value matrix (`arguments.get(0)`). Read through gzip when the name ends with `.gz`, otherwise as plain text.
2. `outputPrefix` — prefix of the output file (`arguments.get(1)`).

**Inputs**

- Value matrix: tab-separated. Lines starting with `#` are skipped; there is no other header handling, so a header line must start with `#`. Columns (0-based index): 0 = chromosome (used only for `-regions`), 1 = start (integer, used only for `-regions`), 2 = end (integer; used for `-regions` and as the row label in the output), 3 to 5 = ignored, 6 and onwards = one numeric value per cell. The text `NA` (any case) is a missing value; every other entry must parse as a number. No BAM is read. Cells are identified only by column position. All rows should have the same number of columns (the code indexes row i with the column count of row j).
- All kept rows go into one matrix regardless of chromosome.
- `-regions` file (optional, `.gz` when the name ends with `.gz`): tab-separated BED; lines starting with `#` are skipped; no first-row skip option. Columns (0-based): 0 chr, 1 start, 2 end, used as written. A matrix row is kept only if a region on the same chromosome (exact name match) overlaps its start..end.
- There are no methylation-pattern options; the tool does not look at reads.

**Options**

| Option | Default | What it does in the code |
|---|---|---|
| `-aveWindow` | 1 | Parsed only. The field is never read; it has no effect. |
| `-permutation` | 0 | Parsed only. The field is never read; it has no effect. |
| `-regions` | null | Keeps only matrix rows that overlap a region in this BED file (same chromosome). |
| `-scoreMode` | 0 | 0: Pearson correlation. 1: variance ratio (the usage text calls it "PoissonBinomial"). 2: `JavaMI.MutualInformation.calculateMutualInformation`. Any other value throws `IllegalArgumentException("no such a mode")` at the first row pair that has more than two shared cells. |
| `-h` | false | See "Help and argument checking": the flag itself is never effective; usage is printed when fewer than 2 tokens are given. |

**Output**

One file: `<outputPrefix>.matrix.txt.gz`, gzip-compressed, tab-separated, UTF-8, no header line. There is one line per kept input row, in input order. For N kept rows each line has 1 + N fields and no trailing tab.

| # | Label | Meaning | Missing value |
|---|---|---|---|
| 1 | (none) | End coordinate (input column 3) of row i, as an integer. The chromosome is not written. | never |
| 2 .. N+1 | (none) | Score of row i against row j, for j = 1..N in input order, including j = i. | `NA` when at most two cells have a value in both rows. `NaN` when the score itself is undefined (for example a Pearson correlation with a constant row). |

**Method**

For rows i and j, the cells with a non-missing value in both rows are selected. If there are two or fewer, the field is `NA`. Otherwise:

- `-scoreMode 0`: commons-math `PearsonsCorrelation.correlation` of the two selected vectors.
- `-scoreMode 1`: `obs / exp`, where `obs` is the commons-math `Variance` of `(value_i + value_j)/2` over the selected cells and `exp` is the mean of `Variance(full row i)` and `Variance(full row j)`. The full rows are passed as stored, including the NaN entries that stand for `NA`. `PoissonBinomial2` is not used. Not verified: commons-math `Variance` with its default constructor is the bias-corrected sample variance, and it returns NaN when the array contains NaN; if so, mode 1 gives `NaN` for every pair that involves a row with a missing value.
- `-scoreMode 2`: `JavaMI.MutualInformation.calculateMutualInformation` on the two selected vectors. Not verified: how JavaMI discretises the values and which logarithm base it uses (library code, not in this repository's source).

No permutation is performed. A `MersenneTwister` is created but not used.

**Behaviour notes**

- The whole matrix is held in memory and all N x N row pairs are computed, including the diagonal and both triangles.
- Progress is logged and the writer flushed every 1,000 output rows.
- Rows from different chromosomes are compared with each other; the output label is only the end coordinate.

**Pitfalls**

- `-aveWindow` and `-permutation` are accepted but do nothing.
- The usage text for `-scoreMode 1` says "PoissonBinomial"; the code computes a variance ratio.
- Two different missing markers occur in the output: `NA` (too few shared cells) and `NaN` (undefined score).
- The output file name suffix `.matrix.txt.gz` is the same as for `TestIndepPoiBinom` and `TestIndepPoiBinomOnBam`; the three formats differ (label column here, none in `TestIndepPoiBinom`, three label columns in `TestIndepPoiBinomOnBam`).
- An invalid `-scoreMode` is detected only after the output file has been opened and the first label written.

**Minimal example** (invocation of `test/unit/stats/SingleCellTestOutputTest.java`)

```
java -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/*" \
  main.java.edu.mit.compbio.bisulfitehic.utils.SingleCellTest -regions /path/to/regions.bed /path/to/matrix.txt /path/to/out/prefix
```

Test input (tab-separated; 4 rows, 5 cells) and the expected default (`-scoreMode 0`, no `-regions`) result shape:

```
#chr  start  end    a  b  c  cell1 cell2 cell3 cell4 cell5
chr1  99     100    .  .  .  1     0     1     0     1
chr1  199    200    .  .  .  1     0     1     0     NA
chr1  299    300    .  .  .  NA    NA    NA    0     1
chr1  65545  65546  .  .  .  0     1     0     1     0
```

The output has 4 lines of 5 fields; the labels are `100`, `200`, `300`, `65546`; the test expects row 1 to hold the scores 1, 1, NA, -1 and row 3 to be all NA. With `-regions` containing `chr1 150 250`, only the row labelled `200` remains.

---

## `TestIndepPoiBinom`

**Purpose** — From the same kind of text matrix (one row per CpG, one value column per cell or allele), it writes a square matrix that compares, for every pair of CpGs or of position windows, the spread across cells of the pairwise-averaged values with the spread of the two units taken separately; with `-permutation` it writes a permutation score instead. Despite the class name, the code computes standard deviations only and does not use `PoissonBinomial2`.

**Command**

```
java -Xmx20G -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/ccInference-0.19-slim.jar:/path/to/bisulfiteHic/lib/gatk-package-distribution-3.3.jar" \
  main.java.edu.mit.compbio.bisulfitehic.utils.TestIndepPoiBinom [options] input.bed.gz outputPrefix
```

Positional arguments, in order:

1. `input.bed.gz` — the value matrix (`arguments.get(0)`). Read through gzip when the name ends with `.gz`, otherwise as plain text.
2. `outputPrefix` — prefix of the output file (`arguments.get(1)`).

**Inputs**

- Value matrix: tab-separated. Lines starting with `#` are skipped; a header line must start with `#`. Columns (0-based index): 2 = position (integer; the third column), 6 and onwards = one numeric value per cell or allele, with `NA` (any case) for missing. Columns 0, 1, 3, 4 and 5 are not read. The chromosome column is ignored, so the file must contain one chromosome only; the code does not check this.
- With `-aveWindow` greater than 1 the rows must be in ascending order of the position column; the code does not check this.
- No BAM is read. Cells or alleles are identified only by column position. There are no methylation-pattern options.

**Options**

| Option | Default | What it does in the code |
|---|---|---|
| `-maxDistance` | -1 | When greater than 0, a pair of units whose positions differ by more than this is written as `NA` without being computed. Values of 0 or less disable the check. |
| `-aveWindow` | 1 | 1: every input row is one unit. Greater than 1: rows are grouped into windows of this many bp on the position column, and the unit is the window (per-column means within the window). |
| `-permutation` | 0 | 0: the field is the ratio `obsStd/expStd`. Greater than 0: that many shuffles are done per pair and the field is `-log10` of the permutation fraction (see Method). |
| `-h` | false | See "Help and argument checking": the flag itself is never effective; usage is printed when fewer than 2 tokens are given. |

**Output**

One file: `<outputPrefix>.matrix.txt.gz`, gzip-compressed, tab-separated, UTF-8, no header line, no row labels and no column labels. For N units there are N lines, each with N fields; every field is followed by a tab, so each line ends with a tab.

| # | Label | Meaning | Missing value |
|---|---|---|---|
| 1 .. N | (none) | Field j of line i is the score of unit i against unit j (units in input order; with `-aveWindow` > 1, windows in ascending order). The diagonal and both triangles are written. | `NA` when the pair is farther apart than `-maxDistance`, or when two or fewer columns have a value in both units. `NaN` or `Infinity` can appear when the ratio is undefined or its denominator is 0. |

The positions of the units are not written. With `-aveWindow` 1, line k and field k correspond to the k-th non-`#` input row. With `-aveWindow` W > 1, they correspond to the k-th non-empty window; windows are `(m*W, (m+1)*W]` on the position column for m = 0, 1, 2, ..., and windows without any input row are left out.

**Method**

- Units with `-aveWindow` 1: each row. Missing cells are stored as NaN. `values_i[z]` is the value of row i in column z.
- Units with `-aveWindow` W > 1: a window starting at `pos` (a multiple of W, beginning at 0) collects the rows with `pos < position <= pos + W`. For each column the window keeps the non-missing values of its rows; `values_i[z]` is their mean (NaN if there are none). The last window of the file is stored at end of input.
- For units i and j, the columns z with a value in both are selected and `a[z] = (values_i[z] + values_j[z]) / 2` is formed. If two or fewer columns are selected the field is `NA`.
- `obsStd` is the commons-math `StandardDeviation` of `a`.
- `expStd` is `(sd_i + sd_j) / 2`, where `sd_i` is `DescriptiveStatistics.getStandardDeviation()` of unit i: with `-aveWindow` 1, of the row's values including the NaN entries; with `-aveWindow` > 1, of all non-missing values of all rows and columns in the window pooled together.
- Without `-permutation`: the field is `obsStd / expStd`.
- With `-permutation` P > 0: the non-missing values (or window column means) of unit i and of unit j are put in two lists; P times, both lists are shuffled with `Collections.shuffle` and the standard deviation of `(i[z] + j[z]) / 2` over the first `min(size_i, size_j)` positions is recorded. The recorded values are sorted, `r = findPosition(sorted, obsStd)` is taken, and the field is `-log10((P - r) / P)`. Because `r` is at most `P - 1`, the largest possible value is `log10(P)`. `expStd` is not used in this mode.
- `Collections.shuffle` is called without a random source, and there is no seed option. A `MersenneTwister` is created but not used.

Not verified: commons-math `StandardDeviation` and `DescriptiveStatistics.getStandardDeviation()` are the bias-corrected sample standard deviation, and return NaN when a value is NaN. If so, with `-aveWindow` 1 and without `-permutation`, every pair that involves a row with an `NA` cell gets `NaN`.

**Behaviour notes**

- All units are held in memory and all N x N pairs are computed.
- With `-aveWindow` 1 the per-column accumulator is still filled for the whole file although it is never used in that mode, which adds memory proportional to the number of values.
- When the permutation fraction is 1 the field is `-log10(1)`, which Java prints as `-0.0`.
- Progress is logged and the writer flushed every 1,000 output lines.

**Pitfalls**

- The output carries no coordinates; the mapping from line number to CpG or window has to be rebuilt from the input as described under Output.
- The chromosome column is ignored; a multi-chromosome file is processed as one coordinate axis.
- `NA`, `NaN` and `Infinity` can all occur as field values.
- Each line ends with a tab; a parser that keeps trailing empty fields sees N + 1 fields.
- The permutation values are not reproducible between runs.
- The output file name suffix is shared with `SingleCellTest` and `TestIndepPoiBinomOnBam`.

**Minimal example** (invocation of `test/unit/stats/TestIndepPoiBinomWindowTest.java`)

```
java -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/*" \
  main.java.edu.mit.compbio.bisulfitehic.utils.TestIndepPoiBinom -aveWindow 10 /path/to/matrix.txt /path/to/out/prefix
```

The test input has six rows (third column 1, 2, 11, 12, 21, 22) and four value columns. Without options the output has 6 lines of 6 fields; with `-aveWindow 10` it has 3 lines of 3 fields (windows (0,10], (10,20], (20,30]).

---

## `TestIndepPoiBinomOnBam`

**Purpose** — Starting from one alignment file per single cell and a VCF of heterozygous SNPs, it assigns read pairs to cell-and-allele groups, bins both mates into fixed position windows, and writes a window-by-window matrix that compares the across-group variance of the methylation level of read pairs linking two windows with the variance within each window; with `-permutation` it writes a permutation score instead.

**Command**

```
java -Xmx20G -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/ccInference-0.19-slim.jar:/path/to/bisulfiteHic/lib/gatk-package-distribution-3.3.jar" \
  main.java.edu.mit.compbio.bisulfitehic.utils.TestIndepPoiBinomOnBam [options] snp_vcf.gz outputPrefix -bams cell1.bam -bams cell2.bam ...
```

Positional arguments, in order:

1. `snp_vcf.gz` — VCF with the SNP genotypes (`arguments.get(0)`).
2. `outputPrefix` — prefix of the output file (`arguments.get(1)`).

`-bams` is formally an option but is required: without it a `CmdLineException` is thrown outside the usage handler.

**Inputs**

- `snp_vcf.gz`: the code requires the tabix index file that htsjdk derives from the VCF name to exist and be readable, otherwise it throws an exception whose message asks for `tabix`. The sequence dictionary is read from the VCF header and passed to a GATK `GenomeLocParser` (which is not used afterwards); the test VCF has contig lines in its header. Not verified: whether a VCF without contig header lines is accepted. All records of all contigs are read into memory. A record is used when it is not filtered, not an indel, not mixed, a point event, and the genotype of the first sample (index 0) is not filtered and heterozygous. Phasing is not checked; the unit test uses an unphased het genotype. Allele A is the first base of the genotype's first allele and allele B the first base of its second allele, in the order of the GT field.
- `-bams`: one file per single cell, given by repeating the option (`-bams a.bam -bams b.bam`). Read groups are not used; the cell is the file. The unit test passes SAM text files. Each file is read sequentially from start to end, two records at a time, and the two records are taken to be the mates of one pair. The file must therefore have the two mates of every pair next to each other and no other records in between (no index is used and coordinate sorting is not suitable). If the two read names differ (ignoring case) a `NumberFormatException` with the message "readname 1 is not equal to readname 2" is thrown. Reads need an `MD` tag (see the common section).
- Chromosome naming: the read's contig name must equal a contig name in the VCF exactly, otherwise the pair is skipped.
- Methylation pattern: fixed CpG; no pattern options.

**Options**

| Option | Default | What it does in the code |
|---|---|---|
| `-bams` | null | Input alignment file of one cell; repeat once per cell. Required. Cell k (0-based order on the command line) gets allele slots `2k` (allele A) and `2k+1` (allele B). |
| `-minBaseQ` | 5 | Minimum base quality of the read base at the SNP and of the read base at each CpG. |
| `-minCellNum` | 2 | A window pair is scored only when more than this many cell-and-allele slots have at least one linking read pair with CpG calls. The comparison is "less than or equal gives NA", and the count is of slots (two per cell), not of cells. |
| `-na` | NA | Text written for a missing value. |
| `-useBadMate` | false | When set, reads without the proper-pair flag are accepted. It takes effect in this class (both in its flag filter and in the methylation call). |
| `-minMapQ` | 30 | Both mates must have at least this mapping quality. |
| `-turnOffBisulfiteFilter` | false | Skips the call to `CcInferenceUtils.bisulfiteIncompleteReads`. |
| `-aveWindow` | 100000 | Window size in bp. A mate is assigned to window `alignment start / aveWindow` (integer division). Windows are always used in this class. |
| `-permutation` | 0 | 0: the field is the ratio `obs/exp`. Greater than 0: that many shuffles are done per window pair and the field is `-log10` of the permutation fraction. |
| `-chr` | chr19 | Used only as the text of the first output column. The chromosome filters that would use it are commented out; reads and SNPs of all chromosomes are processed. |
| `-h` | false | See "Help and argument checking": the flag itself is never effective; usage is printed when fewer than 2 tokens are given. |

**Output**

One file: `<outputPrefix>.matrix.txt.gz`, gzip-compressed, tab-separated, UTF-8, no header line. There is one line per window index from 0 up to the highest window index that received a mate of a SNP-overlapping same-contig pair (W windows in total), including windows without data. Each line has 3 + W fields and no trailing tab.

| # | Label | Meaning | Missing value |
|---|---|---|---|
| 1 | (none) | The value of `-chr`, as text. | never |
| 2 | (none) | Window start, `i * aveWindow`. | never |
| 3 | (none) | Window end, `(i + 1) * aveWindow`. | never |
| 4 .. 3+W | (none) | Score of window i against window j, for j = 0 .. W-1, including j = i. | The `-na` text when too few slots have linking pairs (see `-minCellNum`) or when the score is NaN. `Infinity` is written as is. |

**Method**

1. SNP table: per contig, the kept het SNPs with their two alleles.
2. For each cell file and each record pair on the same contig: the pair is kept if read 1, or otherwise read 2, overlaps a kept SNP; the SNP used is the lowest-coordinate one overlapping read 1, or else read 2. Both mates must pass the flag filter (mapped, primary, mapping quality at least `-minMapQ`, not QC-failed, not duplicate, paired, proper pair unless `-useBadMate`).
3. Strand rule at the SNP: when the two alleles are C and T, pairs whose SNP-carrying read is a C-to-T read (first-of-pair forward or second-of-pair reverse) are skipped; when the two alleles are G and A, pairs whose SNP-carrying read is a G-to-A read are skipped.
4. The read base at the SNP must exist in the read (not deleted or clipped) and have base quality at least `-minBaseQ`. The pair is assigned to allele A if the base matches allele A "in bisulfite space" (`BisulfiteHicUtils.baseMatchInBisulfiteSpace`: equal bases match; additionally reference C matches read T on C-to-T reads and reference G matches read A on G-to-A reads), otherwise to allele B if it matches allele B, otherwise the pair is dropped.
5. CpG calls (0/1) of each mate are taken from `readsMethySummaryAtEachLoc`. A pair in which neither mate has a call is dropped. Mate 1's list of calls is stored under the read name in window `start1 / aveWindow`, slot (cell, allele); mate 2's list in window `start2 / aveWindow`, same slot.
6. For windows i and j and each slot z: over the read names present in slot z of both windows, the calls in window i are summed (`sum_i`, `count_i`) and those in window j likewise. The slot contributes the pooled level `(sum_i + sum_j) / (count_i + count_j)`, and the per-window levels `sum_i / count_i` and `sum_j / count_j`.
7. If the number of contributing slots is less than or equal to `-minCellNum`, the field is the `-na` text.
8. `obs` is the commons-math `Variance` of the pooled levels over slots; `exp` is the mean of the `Variance` of the window-i levels and the `Variance` of the window-j levels. Without `-permutation` the field is `obs / exp` (the `-na` text if NaN).
9. With `-permutation` P > 0: the per-slot `(sum, count)` pairs of window i and of window j are shuffled independently P times with `Collections.shuffle`; each time the variance of `(sum_i + sum_j) / (count_i + count_j)` over the re-paired slots is recorded. With `r = findPosition(sorted recorded values, obs)`, the field is `-log10((P - r) / P)` (the `-na` text if NaN). There is no seed option; the unit test notes that the permutations are not seeded.

The unit test's worked case: two cells, one read pair per allele, pooled levels 0, 0.25, 0.25, 0.5; the expected window 0 / window 2 field is 0.5.

Not verified: commons-math `Variance` with its default constructor is the bias-corrected sample variance (the test's expected value 0.5 = (1/24)/(1/12) is consistent with that).

**Behaviour notes**

- Windows are indexed by position only. Reads from every contig that has SNPs in the VCF fall into the same window list, and the first output column is just the `-chr` text.
- When both mates start in the same window, mate 2's calls replace mate 1's for that read name in that window (both are stored under the same read name in the same table). Diagonal fields are computed from these tables, with the same reads counted on both sides.
- Windows are created for a SNP-overlapping pair before the flag filter is applied, so the matrix size is set by pair positions, not by pairs that pass the filters.
- Only pairs overlapping a het SNP contribute; a pair uses a single SNP.
- No distance filter is applied between mates.
- Each window holds `2 * (number of cells) + 1` slot tables; the last one is never filled.
- All kept read names and their calls are held in memory for all cells; the matrix is W x W.
- After reading, the code logs `methyStatus.get(0).size()`; if no same-contig pair overlapped any het SNP the window list is empty and this call throws `IndexOutOfBoundsException`.
- Progress is logged and the writer flushed every 1,000 output lines.

**Pitfalls**

- `-chr` does not restrict anything. To analyse one chromosome, the input files themselves have to contain only that chromosome's pairs (or only that chromosome's SNPs in the VCF, since pairs on contigs without SNPs are skipped).
- The input alignment files must have mates adjacent. Any extra record (secondary, supplementary, unpaired) shifts the pairing and stops the run with the read-name exception. The code calls `equalsIgnoreCase` on the contig name of the first record of each pair without a null check. Not verified: an unmapped record without a reference name would cause a `NullPointerException` there.
- `-minCellNum 2` means at least three cell-and-allele slots, not two cells.
- The usage text of `-aveWindow` says the default does not use windows; in this class windows are always used and the default is 100000.
- The source comment says only phased het SNPs are considered; the code does not test phasing, and "allele A/B" follow the GT field order of each record.
- The permutation values are not reproducible between runs.
- The output file name suffix is shared with `SingleCellTest` and `TestIndepPoiBinom`.

**Minimal example** (invocation of `test/unit/stats/TestIndepPoiBinomOnBamTest.java`)

```
java -cp "/path/to/bisulfiteHic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfiteHic/lib/*" \
  main.java.edu.mit.compbio.bisulfitehic.utils.TestIndepPoiBinomOnBam -turnOffBisulfiteFilter -aveWindow 10000 -chr chr1 -bams /path/to/cell0.sam -bams /path/to/cell1.sam /path/to/snp.vcf.gz /path/to/out/prefix
```

In the test (mates at about 1 kb and 25 kb, one het SNP), the output has 3 lines of 6 fields; the first line starts with `chr1`, `0`, `10000`, and its sixth field (window 0 against window 2) is 0.5.

---

## Helper classes

## `PoissonBinomial2`

`PoissonBinomial2` is a library class, not a tool. Its `main()` exists only inside a block comment, so the class has no entry point; running `java ... main.java.edu.mit.compbio.bisulfitehic.utils.PoissonBinomial2` cannot start a program. Java answers `Error: Main method not found in class ...PoissonBinomial2`. The commented-out method would have drawn ten random probabilities and printed the PMF and CDF at 10 to standard error; as written it calls `r.nextdouble()`, which is not a `java.util.Random` method.

What the class provides: `new PoissonBinomial2(double[] p)` computes the distribution of the number of successes among independent trials with success probabilities `p` (a recurrence on the odds `p/(1-p)`; trials with `p == 1.0` are removed first and the result is shifted by their count). `get_pmf(k)` returns P(X = k) and `get_cdf(k)` returns P(X <= k); both return NaN when `k` is larger than the number of trials, and neither checks for negative `k`. The only caller in the repository is `MethyCorSepAlleleSC`. None of the five tools in this document uses it, including the two whose names contain "PoiBinom".

## `FisherExactTest`

`FisherExactTest` is a helper class without a `main()` and without options. `new FisherExactTest(maxSize)` precomputes log-factorials up to `maxSize`. For a 2x2 table with cells `a, b, c, d` it offers: `getP` (the hypergeometric probability of exactly this table), `getRightTailedP` (sum over this table and those obtained by increasing `a` and `d`), `getLeftTailedP` (the same in the other direction), `getCumlativeP` (the right tail when `a*d >= b*c`, otherwise the left tail), and `getTwoTailedP` (sum of the probabilities of this table and of every table with the same margins, in either direction, whose probability is at most the observed probability times `1 + 1e-7`). Every method returns NaN when `a+b+c+d` exceeds `maxSize`. The only caller in the repository is `MethyCorAcrossHiccups`, which constructs it with `maxSize` 10000. The unit test checks the two-sided value against exact integer arithmetic for all tables with n up to 30.

---

## Known limitations

Found while reading the source for this guide (November 2026). They are described here, not fixed.

All five tools

1. `-h` does nothing; usage is printed only when too few tokens are given.
2. No permutation can be seeded, so permutation p-values differ between runs.
3. The permutation rank never reaches n, so the p-value never reaches 1 and `-log10` scores are capped at `log10(permutations)`.
4. Missing values are written in several ways: `NA`, `NaN`, possibly `Infinity`, and `-0.0` for a zero permutation score.
5. The three text-matrix tools share the suffix `.matrix.txt.gz` with three different layouts: `TestIndepPoiBinom` writes no labels and a trailing tab on every line, `SingleCellTest` writes one label column (the end coordinate), `TestIndepPoiBinomOnBam` writes chromosome, start and end.

`SingleBaseResHicMatrix`, `SingleBaseResHicMatrixV2`

6. **The CpG list is not a standard BED.** The start column must be the 1-based position of the C (the unit test writes `C_position`, `C_position + 1`). A standard 0-based BED gives only the header line.
7. `-useBadMate` has no effect: reads without the proper-pair flag are always dropped. In BAMs aligned with `-pbat` most long-range pairs lack that flag after `samtools fixmate` (81% in one production cell; see [aligners.md](aligners.md#known-limitations)), so these tools see few Hi-C contacts in such BAMs.
8. `-bgBam` has no effect in either class. In V1 the file is still opened, so a wrong path fails.
9. Both classes write `<prefix>.singleBaseResHicMatrix.bed.gz`; V2 does not use a different name.
10. `-minReadsInWindow` is a per-CpG (V1 only) and shared-read-name threshold; these classes have no windows.
11. With `-bedpe`, all CpGs of a named chromosome are still loaded, and the `name` anchors may come from different BEDPE rows when a CpG overlaps several anchors. A `-bedpe` or `-regions` file from which no row survives leaves the run unrestricted without a message.
12. `-skipProcessFirstRow` in V1 also drops the first row of the `-regions` file.
13. The strand column of the CpG list is read and not used. Supplementary alignments are not excluded.

`SingleCellTest`

14. `-aveWindow` and `-permutation` are declared and never used.
15. `-scoreMode 1` is labelled "PoissonBinomial" but computes a ratio of variances.

`TestIndepPoiBinom`

16. The chromosome column is ignored. With `-aveWindow 1`, rows that contain `NaN` enter the expected standard deviation.

`TestIndepPoiBinomOnBam`

17. `-chr` only sets the text of output column 1; reads from all contigs are pooled into windows by coordinate. Run it on a BAM that holds one chromosome.
18. `-minCellNum` counts cell-and-allele slots (two per BAM) and requires strictly more than the value.
19. Records are consumed two at a time, so mates must be adjacent with no extra record; mate 2's calls overwrite mate 1's when both start in the same window.
20. A run in which no pair overlaps a heterozygous SNP, or without `-bams`, ends in an exception.
21. Phasing is not required; allele A and B follow the order of the GT field.
