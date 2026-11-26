# Methylation concordance tools (`MethyCor*`)

Part of the [usage guide](USAGE.md). These five Java tools ask whether the methylation state at two loci goes together, on the same read pair (bulk) or across single cells.

| Tool | Input | Methylation context | Choose it when |
|---|---|---|---|
| [`MethyCorAcrossHiccups`](#methycoracrosshiccups) | one BAM + BEDPE anchors | any pattern (GCH, HCG, WCG, CG) | you want per-read-pair methylation at anchor pairs (loops); also has shuffled, short-range and allele-specific modes |
| [`MethyCorAcrossHiccupsSC`](#methycoracrosshiccupssc) | per-cell BAMs and bigWigs + BEDPE + CpG list | from bigWigs | single cells, at anchor pairs |
| [`MethyCorEachChrScBw`](#methycoreachchrscbw) | per-cell bigWigs | from bigWigs | single cells, window-by-window matrix of one chromosome |
| [`MethyCorGlobal`](#methycorglobal) | one BAM | CpG only | window-by-window matrix of one chromosome from read pairs |
| [`MethyCorSepAlleleSC`](#methycorsepallelesc) | per-cell BAMs + SNP VCF | CpG only | allele-resolved, single cells |

- Column numbers in the Output sections are 1-based (as `cut`, `awk` and R count them). Column numbers in the Inputs sections are 0-based where the code indexes them that way, and this is said each time.
- `MethyCorAcrossHiccupsSC` and `MethyCorSepAlleleSC` have no unit test. The other three are covered by `test/unit/concordance/`.
- Two helpers live in `lib/ccInference-0.19-slim.jar` without source in this repository (the conversion filter and the IUPAC matcher); what is said about them was read from the compiled classes.
- Sentences that start with "Not verified:" rest on library behaviour that was not checked.

---

## Shared concepts

### Which tools use which helper

| Tool | Reads BAM | Methylation call comes from | Pattern options |
|---|---|---|---|
| `MethyCorAcrossHiccups` | yes | per-read call, `readsMethySummaryGeneral` | yes (`-methyPatternSearch`, `-methyPatternCytPos`, `-bsConvPatternSearch`, `-bsConvPatternCytPos`) |
| `MethyCorGlobal` | yes | per-read call, `readsMethySummary` | none: fixed to `CG` with a `CH` conversion filter |
| `MethyCorSepAlleleSC` | yes | per-cytosine call, `readsMethySummaryAtEachLoc` | none: fixed to `CG` with a `CH` conversion filter |
| `MethyCorAcrossHiccupsSC` | yes, only to find which sites a read pair connects | bigWig values (`-bws`, `-bwCounts`) | none |
| `MethyCorEachChrScBw` | no | bigWig values (`-bws`) | none |

None of the five classes has a second-pattern option (`-methyPattern2ndSearch` / `-methyPattern2ndCytPos` exist only in `LongRangeAsm`, which is outside this section).

### Methylation pattern options (`MethyCorAcrossHiccups` only)

- `-methyPatternSearch` is a string of IUPAC nucleotide codes that is matched against the **reference** sequence of the read (not against the read bases). The matcher accepts the standard codes: `A C G T`, `R`=A/G, `Y`=C/T, `S`=C/G, `W`=A/T, `K`=G/T, `M`=A/C, `B`=C/G/T, `D`=A/G/T, `H`=A/C/T, `V`=A/C/G, `N`=any. Every position of the pattern must match. A reference byte that is not A/C/G/T (an `N`, or the `-` that htsjdk writes for an inserted base) never matches.
- `-methyPatternCytPos` is the 1-based position, inside the pattern, of the cytosine whose methylation is read. Examples from the usage strings: `CG` with 1, `HCG` with 2, `GCH` with 2. The code does not check that the pattern letter at that position is `C`.
- `-bsConvPatternSearch` and `-bsConvPatternCytPos` describe the context used by the bisulfite-conversion filter (see below), with the same syntax. Defaults are `CH` and 1. The README examples for NOMe-HiC use `WCH` and 2.
- The reference sequence of each read is rebuilt from the read's `MD` tag and CIGAR (`SequenceUtil.makeReferenceFromAlignment`). No reference FASTA is read. htsjdk throws `Cannot create reference from SAMRecord with no MD tag` when the tag is missing, so the BAM must have been passed through `samtools calmd` (the README examples use `*.calmd.bam`).
- Soft-clipped bases are removed from read, quality and reference before matching. Inserted bases stay in the arrays (reference byte `-`), deleted reference bases are not represented.

### How the strand of a read is determined

- Only SAM flags are used: `0x10` (reverse strand), `0x1` (paired) and `0x80` (second of pair). No tag (`XG`, `XR`, ...) is read by any of these classes.
- A read is treated as a **C-to-T read** when it is (first of pair or unpaired) and forward, or second of pair and reverse.
- Every other read (first of pair or unpaired and reverse, or second of pair and forward) is treated as a **G-to-A read**.
- A read without the paired flag is always treated like a first-of-pair read.

### How a base is called methylated or unmethylated

- C-to-T read: the pattern is searched on the forward reference. At the pattern's cytosine position, read base `C` counts as methylated, read base `T` counts as unmethylated.
- G-to-A read: in `readsMethySummaryGeneral` the read and the reference are reverse-complemented and the base qualities reversed, then the same rule is applied. On the forward strand this is equivalent to: read `G` at a reference `G` = methylated, read `A` = unmethylated, with the pattern matched on the reverse-complemented reference.
- Any other read base, or a base with quality below `-minBaseQ`, is not counted.
- The fixed-`CG` helpers used by `MethyCorGlobal` and `MethyCorSepAlleleSC` do not reverse-complement. For a C-to-T read they look at reference `C` followed by `G` and count read `C`/`T`. For a G-to-A read they look at reference `G` preceded by `C` and count read `G`/`A`.
- The per-read result is three integers: number of methylated calls (`numC`), number of unmethylated calls (`numT`), and number of pattern sites with a call (`numCg`). `numCg` always equals `numC + numT`. The per-read methylation level written by the tools is `numC / (numC + numT)`.

### Bisulfite-conversion filter

- The filter is on by default and is switched off with `-turnOffBisulfiteFilter` (where that option exists).
- It is `CcInferenceUtils.bisulfiteIncompleteReads`, called with fixed arguments `(…, false, true, true, 0.1, 0.4, 0, pattern, cytPos)`. From the compiled library (its source is not in this repository):
  - For a G-to-A read, read and reference are reverse-complemented first.
  - It counts `U` = number of positions where the **read's own sequence** matches the whole conversion pattern (for `WCH`: the read still shows `C` between a W and an H), and `R` = number of positions where the **reference** matches the whole pattern.
  - The read is rejected when `U / R >= 0.4`. A read with `R = 0` and `U = 0` passes. A read with `R = 0` and `U > 0` is rejected.
  - The mismatch-fraction test (the `0.1` argument) is disabled by the `false` argument.
  - The 5'-conversion-start feature is disabled by the `0` argument, so a read that passes is used over its whole aligned length.
- A rejected read returns zero counts, which the tools then treat as "no methylation call".
- In `MethyCorGlobal` and `MethyCorSepAlleleSC` the filter pattern is hard-coded to `CH`, position 1.

### Default read filters

A record is dropped when any of the following is true. The same tests are applied twice: when reads are collected from a region (`failFlagFilter` in each tool) and again inside the per-read methylation helper.

| Test | Flag / field | Default | Option that changes it |
|---|---|---|---|
| unmapped | `0x4` | dropped | none |
| secondary alignment | `0x100` | dropped | none |
| duplicate | `0x400` | dropped | none |
| fails vendor QC | `0x200` | dropped | none |
| mapping quality | MAPQ `<` threshold | 30 in `MethyCorAcrossHiccups`, `MethyCorGlobal`; 20 in `MethyCorAcrossHiccupsSC`, `MethyCorSepAlleleSC` | `-minMapQ` |
| not paired | `0x1` not set | dropped | `-useUnpaired` (only `MethyCorAcrossHiccups`, `MethyCorGlobal`) |
| "bad mate" | proper-pair flag `0x2` not set | dropped | `-useBadMate` |
| base quality (per base, in the methylation call and at the SNP base) | QUAL `<` threshold | 5 | `-minBaseQ` (not in `MethyCorAcrossHiccupsSC`, which makes no base-level call) |

- `-useBadMate` means "also keep records whose proper-pair flag is not set". All README Analysis 4 commands pass it. In the unit test (`global` scenario) long-range pairs flagged not-proper give no result without it.
- The supplementary flag (`0x800`) is not tested anywhere in these classes.
- Duplicates are recognised only by the flag, so the BAM must be duplicate-marked beforehand.
- `-useUnpaired` must be given together with `-useBadMate`: with `-useUnpaired` alone the code calls `getProperPairFlag()` on an unpaired record, and htsjdk throws `IllegalStateException: Inappropriate call if not paired read`.

### How the two ends of a pair are found

- The tools query the BAM by region and match the two ends **by read name** in a `HashMap`. The two mates must therefore have identical read names (or names that become identical after `-removeReadnameSuffix` cuts the last two characters).
- In most modes a name seen twice inside the same region is removed again ("both ends in the same anchor"). The logic is a toggle: a third record with the same name would add the name back.
- "End 1" and "end 2" in the outputs mean "the record found in region 1" and "the record found in region 2". They do not mean read 1 and read 2 of the sequencer.

### How the distance between the two ends is measured

- `BisulfiteHicUtils.mateDistance(r)` returns `|POS − PNEXT|` (alignment start of the record minus alignment start of its mate, absolute value) when both are on the same contig.
- It returns `Integer.MAX_VALUE` when the read is unpaired, the mate is unmapped, or the mate is on another contig.
- TLEN is not used.
- This distance is used only by `-onlyShort` / `-shortLowerBound` in `MethyCorAcrossHiccups` and by `-onlyShort` in `MethyCorAcrossHiccupsSC`. In the default long-range mode no read-level distance is tested: only the distance between the two BEDPE anchors is tested (`-minDistEnds`).

### Common command-line behaviour

- Main classes are `main.java.edu.mit.compbio.bisulfitehic.utils.<ClassName>`.
- Options and positional arguments may be mixed in any order (the README puts positionals first, the unit test puts options first).
- List options (`-bams`, `-bws`, `-bwCounts`) take one value per occurrence and are repeated: `-bams a.bam -bams b.bam`.
- Running a class with too few command-line tokens prints the usage string and the option list to stderr and returns.
- Progress messages are written through log4j `BasicConfigurator` (console).
- Jars, verified by listing jar contents:
  - `target/bisulfitehic-0.38-jar-with-dependencies.jar`: the tool classes, htsjdk, args4j, commons-math3, commons-lang3, log4j.
  - `lib/ccInference-0.19-slim.jar`: `CcInferenceUtils`, `BaseUtilsMore`.
  - `lib/gatk-package-distribution-3.3.jar`: `GenomeLocParser`, `BaseUtils`, `org.apache.commons.lang.ArrayUtils`.
  - `lib/java-genomics-io.jar`: `edu.unc.genomics.io.BigWigFileReader` (needed by the two bigWig tools).
  - `lib/igv.jar`: `org.broad.igv.bbfile.BBFileReader`, which `BigWigFileReader` calls.
- The README gives a classpath (the first three jars) only for `MethyCorAcrossHiccups`. `test/unit/concordance/run.sh` uses the fat jar plus `lib/*` plus `jbwa/jbwa-1.0.0/jbwa.jar`.
- Not verified: that the five-jar classpath shown below for the bigWig tools is sufficient at run time (class locations were checked, the commands were not run).

---

## `MethyCorAcrossHiccups`

**Purpose** — For each row of a BEDPE file (a pair of anchors, for example chromatin loops) it finds read pairs with one end in each anchor and writes, per read pair, the methylation level of each end. The correlation itself is computed afterwards by the user (the README has an R example). Choose it for bulk or pseudobulk BAMs when you want read-level concordance at given anchor pairs; it is the only tool here with configurable methylation patterns (GCH, HCG, WCG, ...). It also has a shuffled-control mode, a short-range mode and an allele-specific mode (`-vcf`).

**Command**

```
java -Xmx20G -cp "<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar:<repo>/lib/ccInference-0.19-slim.jar:<repo>/lib/gatk-package-distribution-3.3.jar" \
  main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorAcrossHiccups [opts] outputPrefix anchors.bedpe[.gz] input.bam
```

Positional arguments:
1. `outputPrefix` — prefix of the output file (`arguments.get(0)`).
2. `hiccups.bedpe.bed[.gz]` — BEDPE file of anchor pairs (`arguments.get(1)`). It is read as gzip only when the name ends in `.gz`.
3. `input.bam` — BAM file (`arguments.get(2)`).

The usage check counts all command-line tokens, including options, and requires at least 3.

**Inputs**

- BAM
  - Queried by region, so it must be coordinate-sorted and indexed. Not verified: the exact error htsjdk gives for a missing index (library behaviour, not read here); the unit test writes a coordinate-sorted BAM with an index.
  - Every record needs an `MD` tag (see Shared concepts).
  - Mates must share the read name. Duplicates must be flagged.
  - Opened with validation stringency `SILENT`.
- BEDPE (tab-separated)
  - Columns used, 0-based: 0 `chr1`, 1 `start1`, 2 `end1`, 3 `chr2`, 4 `start2`, 5 `end2`.
  - Unless `-useGeneralBedPe` is given, columns 7, 11 and 15 (0-based) are parsed as numbers and copied to the output as `ob`, `exp`, `fdr`. The file then needs at least 16 columns with numbers in those three.
  - Lines starting with `#` are skipped. With `-skipProcessFirstRow` the first physical line is skipped too. Any other header line makes `Integer.parseInt` fail.
  - `start`/`end` are passed unchanged to htsjdk's `queryOverlapping` (and to the VCF query). htsjdk treats them as 1-based inclusive (library contract); no +1 is added to the BED start.
  - Chromosome names must match the BAM (and VCF) names exactly; they are passed through as written.
- VCF (only with `-vcf`)
  - bgzip-compressed with a tabix index next to it; the tool stops with `... file's index ... does not exist, please use tabix to index it` otherwise.
  - The sequence dictionary is taken from the VCF header. Not verified: behaviour when the header has no contig lines.
  - Only records that are unfiltered, not indels, not mixed, single-position, and heterozygous with an unfiltered genotype in sample number `-sample` are used. Allele A is the first allele of the genotype, allele B the second (first base of each).

**Options**

| Option | Default | What it does in the code |
|---|---|---|
| `-minBaseQ` | 5 | Minimum base quality for a methylation call and (ASM mode) for the SNP base. |
| `-turnOffBisulfiteFilter` | false | Skips the conversion filter. |
| `-useBadMate` | false | Keeps records without the proper-pair flag. |
| `-useUnpaired` | false | Keeps records without the paired flag. Needs `-useBadMate` as well (see Shared concepts). Has no useful effect with `-onlyShort`, because an unpaired read has mate distance `Integer.MAX_VALUE`. |
| `-removeReadnameSuffix` | false | Cuts the last two characters of every read name before matching ends. Not applied in ASM mode (`-vcf`). |
| `-minMapQ` | 30 | Minimum mapping quality. |
| `-minDistEnds` | 20000 | BEDPE row filter: a same-chromosome row is skipped when `|start1−start2|` or `|end1−end2|` is below this value. Not applied when `-onlyShort` is 0 or larger. It is not a read-level distance. |
| `-vcf` | null | Switches to ASM mode: writes `<prefix>.asm.txt.gz` instead of `<prefix>.methyCor.txt.gz`. Takes precedence over `-randomHiccupsInput` and `-onlyShort`. |
| `-coverageRef` | 2 | ASM mode only: minimum number of usable reads carrying allele A at the SNP. |
| `-coverageAlt` | 2 | ASM mode only: minimum number of usable reads carrying allele B at the SNP. |
| `-minCgNum` | 1 | Minimum number of pattern sites with a call in a read (each end). Applied in every mode, not only ASM, except for rows written with `-outputCount`. |
| `-sample` | 0 | ASM mode only: 0-based index of the VCF sample whose genotype is used. |
| `-randomHiccupsInput` | false | Shuffled-control mode (see Behaviour notes). Takes precedence over `-onlyShort`. |
| `-randomNotLinked` | false | Only with `-randomHiccupsInput`: pair any read of anchor 1 with any read of anchor 2, whether or not their mates are in the other anchor. |
| `-shortLowerBound` | -1 | Only in `-onlyShort` mode: records with mate distance below this value are dropped. |
| `-onlyShort` | -1 | When > 0: short-range mode; records with mate distance above this value are dropped, and the `-minDistEnds` row filter is off. When exactly 0: the `-minDistEnds` row filter is off but the normal long-range mode still runs. |
| `-randomSampleReads` | -1 | With `-randomHiccupsInput` (without `-randomNotLinked`): number of shuffle rounds per BEDPE row (1 round when not set). With `-onlyShort`: stop after `N+1` rows per BEDPE row. No effect in the default mode, with `-randomNotLinked`, or in ASM mode. It never draws a random subset. |
| `-minReads` | 0 | Default and short mode: a BEDPE row is skipped when it has fewer than this many name-matched pairs. Random mode: the row is used only when the count is strictly greater than this value. Not used in ASM mode. |
| `-skipProcessFirstRow` | false | Skips the first physical line of the BEDPE file. |
| `-diffChr` | false | Keeps rows whose two anchors are on different chromosomes. In `-onlyShort` mode the second chromosome is ignored. |
| `-useGeneralBedPe` | false | Does not read columns 7/11/15; `ob`, `exp`, `fdr` are written as `NaN`. |
| `-outputCount` | false | Default long-range mode only: writes one row for every name-matched pair, including pairs without any methylation call, and appends a 16th column with the total number of read pairs in the BAM. Ignored in random, short and ASM modes. |
| `-methyPatternSearch` | `CG` | IUPAC pattern whose cytosine is called. |
| `-methyPatternCytPos` | 1 | 1-based cytosine position in `-methyPatternSearch`. |
| `-bsConvPatternSearch` | `CH` | IUPAC pattern used by the conversion filter. |
| `-bsConvPatternCytPos` | 1 | 1-based cytosine position in `-bsConvPatternSearch`. |
| `-h` | false | No effect of its own. The usage text is printed only when fewer than 3 tokens are given; `-h` together with 3 or more tokens runs the tool. |

**Output**

One gzip file, no header line, tab-separated. Row order inside a BEDPE row follows `HashMap` iteration and is not sorted.

`<outputPrefix>.methyCor.txt.gz` (all modes except `-vcf`): one row per read pair.

| Col | Name | Meaning |
|---|---|---|
| 1 | chr1 | Anchor 1 chromosome, from the BEDPE row. |
| 2 | start1 | Anchor 1 start, from the BEDPE row (not the read position). |
| 3 | end1 | Anchor 1 end, from the BEDPE row. |
| 4 | chr2 | Anchor 2 chromosome (in `-onlyShort` mode: chr1 again). |
| 5 | start2 | Anchor 2 start, from the BEDPE row. |
| 6 | end2 | Anchor 2 end, from the BEDPE row. |
| 7 | region_name | `chr1:start1:end1:chr2:start2:end2`. In `-onlyShort` mode: `chr1:start1:end1:start2:end2` (5 fields). |
| 8 | id | Running number of the row within the BEDPE row. Starts at 1 in the default mode and at 0 in random and short modes. |
| 9 | ob | BEDPE column 7 (0-based), or `NaN` with `-useGeneralBedPe`. |
| 10 | exp | BEDPE column 11 (0-based), or `NaN`. |
| 11 | fdr | BEDPE column 15 (0-based), or `NaN`. |
| 12 | methy_end1 | `numC/(numC+numT)` of the read in anchor 1, fraction 0–1. With `-outputCount`, `NaN` when the read has no call. |
| 13 | count_end1 | `numC+numT` of that read (number of pattern cytosines with a call). |
| 14 | methy_end2 | Same as column 12 for the read in anchor 2. |
| 15 | count_end2 | Same as column 13 for the read in anchor 2. |
| 16 | total_read_pairs | Only with `-outputCount` in the default mode: number of records in the whole BAM that pass the flag filter, divided by 2 (integer division). Same value on every row. |

`<outputPrefix>.asm.txt.gz` (only with `-vcf`): one row per heterozygous SNP per anchor that passes the coverage filters. 28 columns.

| Col | Name | Meaning |
|---|---|---|
| 1 | snp_chr | SNP chromosome. |
| 2 | snp_start | SNP position minus 1 (0-based). |
| 3 | snp_end | SNP position (1-based). |
| 4 | alleleA | First base of the first genotype allele. |
| 5 | alleleB | First base of the second genotype allele. |
| 6 | region_name | `chr1:start1:end1:start2:end2` (5 fields, chr2 is not included). |
| 7 | id | 0-based counter of SNPs written for this anchor of this BEDPE row. |
| 8–10 | ob, exp, fdr | As above. |
| 11 | A_all_methy | Methylated calls summed over all usable reads carrying allele A at the SNP. |
| 12 | A_all_unmethy | Unmethylated calls, same reads. |
| 13 | B_all_methy | Methylated calls, reads carrying allele B. |
| 14 | B_all_unmethy | Unmethylated calls, reads carrying allele B. |
| 15 | diff_all | Methylation of B minus methylation of A, in percentage points (×100). `NaN` when either allele has no call. |
| 16 | p_all | Two-sided p-value of the 2×2 table (columns 11–14). Fisher exact test when the table total is at most 10000, chi-square test above that. |
| 17–20 | A_methy, A_unmethy, B_methy, B_unmethy (SNP end, linked) | Same counts restricted to SNP reads whose read name is also found in the other anchor. |
| 21 | diff_snp_end | As column 15, for columns 17–20. |
| 22 | p_snp_end | As column 16, for columns 17–20. |
| 23–26 | A_methy, A_unmethy, B_methy, B_unmethy (other end) | Calls of the mate found in the other anchor, grouped by the allele of the SNP read. |
| 27 | diff_other_end | As column 15, for columns 23–26. |
| 28 | p_other_end | As column 16, for columns 23–26. |

Missing values are written as `NaN`.

**Behaviour notes**

- Mode selection, in this order: `-vcf` → ASM mode; else `-randomHiccupsInput` → random mode; else `-onlyShort > 0` → short mode; else default long-range mode.
- Row filter before any mode: a row is skipped when the anchors are on different chromosomes (unless `-diffChr`), or when they are on the same chromosome and either `|start1−start2|` or `|end1−end2|` is below `-minDistEnds` (unless `-onlyShort >= 0`).
- Default mode: all records overlapping anchor 1 and all records overlapping anchor 2 are collected (flag filter, name toggle). A pair is every read name present in both collections. The row is written when both ends have at least one call and at least `-minCgNum` sites.
- Random mode without `-randomNotLinked`: the name-matched pairs are collected as in the default mode, then the anchor-2 reads are shuffled among the pairs and each anchor-1 read is written with the anchor-2 read of another (randomly assigned) pair. Nothing prevents a read from being assigned its own mate. `-randomSampleReads N` repeats the shuffle N times, so up to N × (number of pairs) rows are written per BEDPE row.
- Random mode with `-randomNotLinked`: the read lists of the two anchors are truncated to the shorter length, the anchor-1 list is shuffled once and the lists are paired position by position.
- Short mode: both anchors are queried on `chr1`. A pair is formed when both mates overlap the union of the two anchors and both have a mate distance in `[-shortLowerBound, -onlyShort]`. "End 1" is the mate met first, "end 2" the other. With anchors far apart, the pairs reported are therefore those with both mates inside the same anchor. A record is never paired with itself in this mode.
- ASM mode: for every heterozygous SNP inside an anchor, reads overlapping the SNP are assigned to allele A or B in bisulfite space (a `C` allele also matches read `T` on a C-to-T read, a `G` allele also matches read `A` on a G-to-A read). For C/T SNPs only G-to-A reads are used; for A/G SNPs only C-to-T reads are used. A SNP is written when allele A has at least `-coverageRef` reads, allele B at least `-coverageAlt` reads, and each of the three count groups (all, SNP end linked, other end) has at least one call.
- Randomness: shuffling uses `Collections.shuffle` without a seed. There is no seed option, so random-mode output is not reproducible between runs.
- Memory and run time: every record overlapping an anchor is held in memory while the row is processed. `-outputCount` reads the whole BAM once before the first row. The output is flushed every 1000 processed rows.

**Pitfalls**

- Without `-useBadMate`, records lacking the proper-pair flag are dropped; the README passes `-useBadMate` in every example.
- Without `-useGeneralBedPe`, a plain 6-column BEDPE fails when columns 7/11/15 are read.
- Columns 2–3 and 5–6 of the output are the anchor coordinates, not read coordinates; the read positions are not written.
- `-randomSampleReads` does not subsample at random in any mode (see option table).
- In the default and random modes there is no check that the anchor-1 record and the anchor-2 record are different alignments. If the two anchors overlap (possible when `-minDistEnds` is lowered, when `-onlyShort 0` is used, or when anchors are wider than `-minDistEnds`), a read lying in both anchors whose mate is in neither is paired with itself.
- Supplementary alignments are not filtered; a third record with the same read name in an anchor re-adds a name that the toggle had removed.
- `-useUnpaired` without `-useBadMate` ends with an `IllegalStateException` on the first unpaired record.
- `-minReads` is "at least" in default/short mode and "more than" in random mode.

**Minimal example** (README Analysis 4.1, GCH with WCH conversion filter)

```
java -Xmx20G -cp "/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorAcrossHiccups GCH_concordance.sample loops.bedpe sample.calmd.bam -methyPatternSearch GCH -methyPatternCytPos 2 -bsConvPatternSearch WCH -bsConvPatternCytPos 2 -useBadMate -useGeneralBedPe
```

Output: `GCH_concordance.sample.methyCor.txt.gz`. Short-range variant (README 4.3): add `-onlyShort 1000 -randomSampleReads 100`. Shuffled control (README 4.4): add `-randomHiccupsInput -randomSampleReads 100`.

---

## `MethyCorAcrossHiccupsSC`

**Purpose** — Single-cell version of the anchor-pair analysis. For each BEDPE row it uses one BAM per cell to find which pairs of listed sites (for example CpGs) are connected by a read pair in that cell, keeps site pairs connected in at least `-minCellNum` cells, and writes for each such cell the methylation value and coverage of both sites taken from that cell's bigWig files. Choose it when you have per-cell BAMs plus per-cell methylation and coverage bigWigs; the methylation context is decided by the site list and the bigWigs, not by options.

**Command**

```
java -Xmx20G -cp "<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar:<repo>/lib/ccInference-0.19-slim.jar:<repo>/lib/gatk-package-distribution-3.3.jar:<repo>/lib/java-genomics-io.jar:<repo>/lib/igv.jar" \
  main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorAcrossHiccupsSC [opts] outputPrefix anchors.bedpe[.gz] cpg_list.bed[.gz] \
  -bams cell1.bam -bws cell1.methy.bw -bwCounts cell1.cov.bw  -bams cell2.bam -bws cell2.methy.bw -bwCounts cell2.cov.bw ...
```

Positional arguments:
1. `outputPrefix` (`arguments.get(0)`).
2. `hiccups.bedpe.bed[.gz]` — BEDPE of anchor pairs (`arguments.get(1)`); gzip only when the name ends in `.gz`.
3. `cpg_list.bed[.gz]` — BED file of the sites to pair (`arguments.get(2)`); gzip only when the name ends in `.gz`.

`-bams`, `-bws` and `-bwCounts` are required, must be non-empty and must have the same number of entries. The i-th entries of the three lists belong to the same cell. If they are missing or of unequal length a `CmdLineException` is thrown (a stack trace, not the usage text).

**Inputs**

- BAMs (`-bams`): one per cell, coordinate-sorted and indexed (region queries). Mates must share the read name. No `MD` tag is needed, because no base-level call is made from the BAM.
- Methylation bigWigs (`-bws`): one per cell. The value at a site is read as the mean of the bigWig over the single 1-based position `site_end` (the BED end column of the site list).
- Coverage bigWigs (`-bwCounts`): one per cell. The value at a site is read as the sum over the same single position. A site counts as covered in a cell when this value is neither `NaN` nor 0.
- Site list BED (tab-separated): columns 0 `chr`, 1 `start` (0-based), 2 `end`. Lines starting with `#` and lines with fewer than 3 columns are skipped. Column 5 (strand) is parsed and stored but never used. Each site is stored as the interval `start+1 .. end` and identified in the output by its `end`.
- BEDPE: same parsing as `MethyCorAcrossHiccups` (columns 0–5; columns 7/11/15 unless `-useGeneralBedPe`; `#` lines skipped).
- Chromosome names must be identical in the BEDPE, the site list, the BAMs and the bigWigs (the site list is looked up with an exact, case-sensitive match).

**Options**

| Option | Default | What it does in the code |
|---|---|---|
| `-bams` | null | Required, repeatable. One BAM per cell. |
| `-bws` | null | Required, repeatable. Methylation-level bigWig per cell, same order as `-bams`. |
| `-bwCounts` | null | Required, repeatable. Coverage (C+T count) bigWig per cell, same order as `-bams`. |
| `-minCellNum` | 2 | A site pair is written only when at least this many cells have a read pair connecting the two sites. |
| `-useBadMate` | false | Keeps records without the proper-pair flag. |
| `-minMapQ` | 20 | Minimum mapping quality. |
| `-minDistEnds` | 20000 | BEDPE row filter, as in `MethyCorAcrossHiccups` (off when `-onlyShort >= 0`). |
| `-randomHiccupsInput` | false | Control mode: each anchor-1 read is paired with a randomly drawn anchor-2 read of a different name, per cell. Used only when both anchors have more than 10 reads in that cell. Takes precedence over `-onlyShort`. |
| `-onlyShort` | -1 | When > 0: records with mate distance above this value are dropped and the row distance filter is off. When 0: only the row distance filter is switched off. |
| `-skipProcessFirstRow` | false | Skips the first physical line of the BEDPE file. |
| `-useGeneralBedPe` | false | Does not read BEDPE columns 7/11/15; `ob`, `exp`, `fdr` are `NaN`. |
| `-h` | false | No effect of its own (see `MethyCorAcrossHiccups`). |

**Output**

`<outputPrefix>.methyCor.txt.gz` — gzip, tab-separated, no header. One row per (site pair, cell) for every cell in which the pair is connected. Row order is not sorted.

| Col | Name | Meaning |
|---|---|---|
| 1 | chr | Chromosome of site 1 (anchor 1 side). |
| 2 | site1_start | `site1_end − 1`. |
| 3 | site1_end | BED end of site 1 (1-based position used for the bigWig lookup). |
| 4 | chr | Same chromosome again. |
| 5 | site2_start | `site2_end − 1`. |
| 6 | site2_end | BED end of site 2. |
| 7 | region_name | `chr:start1:end1:start2:end2` of the BEDPE row. |
| 8 | cell | 0-based index of the cell in the `-bams` list. |
| 9–11 | ob, exp, fdr | BEDPE columns 7/11/15 (0-based) or `NaN`. |
| 12 | methy1 | Mean of the cell's methylation bigWig at site 1. `NaN` when the bigWig has no value there. |
| 13 | count1 | Sum of the cell's coverage bigWig at site 1. |
| 14 | methy2 | As column 12 for site 2. |
| 15 | count2 | As column 13 for site 2. |

**Behaviour notes**

- Rows with anchors on different chromosomes are always skipped (there is no `-diffChr`).
- Per BEDPE row and per cell: records overlapping each anchor are collected (flag filter, name toggle, unpaired reads always dropped). For each name found in both anchors, all listed sites overlapped by the anchor-1 record are combined with all listed sites overlapped by the anchor-2 record, keeping only sites with non-zero coverage in that cell's coverage bigWig.
- The methylation values come from the bigWigs. They are cell-level values at the site, not the call made on the connecting read.
- A record is not paired with itself in the default and short modes.
- Short mode (`-onlyShort N`): unlike the bulk tool, the two anchors are still queried separately and a name with both mates inside the same anchor is removed. A pair is counted only when one mate overlaps anchor 1 and the other overlaps anchor 2.
- Randomness: the random mode uses an unseeded `MersenneTwister`; there is no seed option.
- Memory and run time: the whole site list is loaded into per-chromosome interval trees. Every BAM is opened and closed again for every BEDPE row. Each candidate site triggers single-base bigWig queries.

**Pitfalls**

- There is no methylation-pattern or base-quality option; what is measured is fixed by the site list and by how the bigWigs were made.
- The bigWig value must sit at the 1-based position equal to the BED end of the site. A site list with 2-bp entries points at the second base.
- `-onlyShort` with identical or far-apart anchors yields nothing (see Behaviour notes).
- A site pair seen in fewer than `-minCellNum` cells (default 2) is dropped, so a single-cell run with one BAM produces no output unless `-minCellNum 1` is given.
- Not verified: what the BigWig reader does when a chromosome is absent from a bigWig file.

**Minimal example** (GCH sites; the GCH context is selected by the site list and the bigWigs, the code has no pattern option)

```
java -Xmx20G -cp "/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar:/path/to/bisulfitehic/lib/java-genomics-io.jar:/path/to/bisulfitehic/lib/igv.jar" main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorAcrossHiccupsSC sc_GCH loops.bedpe GCH_sites.bed.gz -bams cell1.bam -bws cell1.GCH.methy.bw -bwCounts cell1.GCH.cov.bw -bams cell2.bam -bws cell2.GCH.methy.bw -bwCounts cell2.GCH.cov.bw -useBadMate -useGeneralBedPe
```

No README or test invocation exists for this class; the command shape above follows the usage string and the option declarations.

---

## `MethyCorEachChrScBw`

**Purpose** — For one chromosome it cuts the chromosome into fixed windows and writes a window × window matrix of Pearson correlations, where each correlation is computed across the given bigWig files (one value per bigWig = mean signal of that file in the window). It uses no BAM and no read pairs. Choose it to ask whether two windows co-vary in methylation across cells (or across alleles) from bigWig tracks alone.

**Command**

```
java -Xmx20G -cp "<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar:<repo>/lib/ccInference-0.19-slim.jar:<repo>/lib/gatk-package-distribution-3.3.jar:<repo>/lib/java-genomics-io.jar:<repo>/lib/igv.jar" \
  main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorEachChrScBw [opts] outputPrefix genome.chrom.sizes chr -bws a.bw -bws b.bw ...
```

Positional arguments:
1. `outputPrefix` (`arguments.get(0)`).
2. `mm9.chroms.size` — chromosome sizes file (`arguments.get(1)`); any genome, the name in the usage string is only an example.
3. `chr` — the chromosome to process (`arguments.get(2)`).

`-bws` is required (at least one). A missing `-bws` throws a `CmdLineException`.

**Inputs**

- Chromosome sizes file: plain text, tab-separated, column 0 name, column 1 length. Every line whose name equals `chr` ignoring case is processed. It is not read as gzip.
- bigWig files (`-bws`): one per cell (or per cell and allele). The `chr` argument is passed to the bigWig reader exactly as typed.

**Options**

| Option | Default | What it does in the code |
|---|---|---|
| `-bws` | null | Required, repeatable. bigWig files; each contributes one value per window (the mean over bases with data). The usage text mentions `-bwCounts`, which does not exist in this class (it is commented out). |
| `-aveWindow` | 1000000 | Window size in bp. Always used; the usage text "default value do not use window" does not describe the code. |
| `-random` | 0 | When > 0: for each window pair the window-2 values are shuffled across bigWigs this many times and the mean of the resulting correlations is written. Each shuffle starts from the previous shuffled order. |
| `-seed` | -1 | Seed of `java.util.Random` for `-random`. `-1` means unseeded. |
| `-shuffleWithinCell` | false | Swaps the window-2 values of consecutive bigWig pairs (1st↔2nd, 3rd↔4th, ...) before correlating. Needs an even number of `-bws`. When set, `-random` and `-seed` are ignored. |
| `-h` | false | No effect of its own (usage is printed only when fewer than 3 tokens are given). |

**Output**

`<outputPrefix>.methyCor.txt.gz` — gzip, no header, no coordinate columns.

- One line per window (row = window 1), one value per window (column = window 2). Each value is followed by a tab, so every line ends with a tab before the newline.
- Row i / column j (1-based) is the window starting at `(i−1) × aveWindow`.
- Windows start at 0, `aveWindow`, 2×`aveWindow`, ... as long as `start < chromosome_length − aveWindow`. The last window of the chromosome (partial, or full when the length is an exact multiple) is not included.
- Value: Pearson correlation coefficient between the vector of window-1 means and the vector of window-2 means over all bigWigs. Missing or undefined values are written as `NaN`.
- The full square matrix is written, including the diagonal and both triangles.

**Behaviour notes**

- A window mean is `queryStats(chr, start, end).getMean()`: the mean over bases that have a bigWig value, with `start` and `end` passed as given (the reader treats the interval as 1-based inclusive, so neighbouring windows share one boundary base).
- bigWigs without data in a window give `NaN` for that window; such entries are not removed before the correlation.
- Not verified: that a `NaN` entry makes the whole correlation `NaN` (expected from commons-math `SimpleRegression`, not run).
- commons-math throws when fewer than 2 bigWigs are given (array length < 2).
- Randomness: only with `-random`; reproducible when `-seed` is set. The generator is created once, so the result of a window pair depends on all pairs processed before it.
- Run time: every window pair queries every bigWig twice; nothing is cached. The number of queries is `windows² × 2 × number of bigWigs`. The output is flushed after each row.

**Pitfalls**

- If `chr` is not in the sizes file, or the chromosome is not longer than `-aveWindow`, the output file is created but empty.
- The output has no labels; window coordinates must be reconstructed from `-aveWindow`.
- `-seed -1` cannot be used as an actual seed.
- Not verified: how the BigWig reader handles the first window, which is queried with start 0.

**Minimal example**

```
java -Xmx20G -cp "/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar:/path/to/bisulfitehic/lib/java-genomics-io.jar:/path/to/bisulfitehic/lib/igv.jar" main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorEachChrScBw GCH_chr1_1Mb genome.chrom.sizes chr1 -bws cell1.GCH.methy.bw -bws cell2.GCH.methy.bw -bws cell3.GCH.methy.bw -aveWindow 1000000
```

The context (GCH here) is whatever the bigWigs contain. No README invocation exists; the unit test only checks the odd-`-bws` rejection of `-shuffleWithinCell`.

---

## `MethyCorGlobal`

**Purpose** — For one chromosome and one BAM it cuts the chromosome into fixed windows and writes a window × window matrix. Each entry is the Pearson correlation, over all read pairs with one end in window 1 and the other in window 2, between the CpG methylation level of the two ends. Choose it for an unbiased chromosome-wide map when no anchor list is available. It measures `CG` only.

**Command**

```
java -Xmx20G -cp "<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar:<repo>/lib/ccInference-0.19-slim.jar:<repo>/lib/gatk-package-distribution-3.3.jar" \
  main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorGlobal [opts] outputPrefix genome.chrom.sizes chr input.bam
```

Positional arguments:
1. `outputPrefix` (`arguments.get(0)`).
2. `mm9.chroms.size` — chromosome sizes file (`arguments.get(1)`).
3. `chr` — chromosome to process (`arguments.get(2)`).
4. `input.bam` (`arguments.get(3)`).

The usage check requires at least 4 command-line tokens.

**Inputs**

- BAM: coordinate-sorted, indexed, with `MD` tags, duplicate-marked, mates sharing the read name (same requirements as `MethyCorAcrossHiccups`).
- Chromosome sizes file: plain text, tab-separated, column 0 name (compared with `chr` ignoring case), column 1 length. `chr` is passed to the BAM query exactly as typed.

**Options**

| Option | Default | What it does in the code |
|---|---|---|
| `-minBaseQ` | 5 | Minimum base quality for a CpG call. |
| `-turnOffBisulfiteFilter` | false | Skips the conversion filter (which is fixed to pattern `CH`, position 1). |
| `-useBadMate` | false | Keeps records without the proper-pair flag. |
| `-useUnpaired` | false | Keeps records without the paired flag. Needs `-useBadMate` as well. |
| `-removeReadnameSuffix` | false | Cuts the last two characters of every read name before matching ends. |
| `-minMapQ` | 30 | Minimum mapping quality. |
| `-minCgNum` | 1 | Minimum number of CpGs with a call in each end. The usage text says "For ASM mode"; this class has no ASM mode. |
| `-randomHiccupsInput` | false | Control mode: each window-1 read is paired with a randomly drawn window-2 read of a different name. Used only when both windows hold more than 10 reads; otherwise the entry is `NaN`. |
| `-randomSampleReads` | -1 | Only with `-randomHiccupsInput`: stop after this many usable random pairs per window pair. No effect otherwise. |
| `-aveWindow` | 1000000 | Window size in bp. |
| `-h` | false | No effect of its own (usage is printed only when fewer than 4 tokens are given). |

**Output**

`<outputPrefix>.methyCor.txt.gz` — gzip, no header, no coordinate columns.

- Same layout as `MethyCorEachChrScBw`: one line per window 1, one tab-terminated value per window 2, full square matrix, windows starting at multiples of `-aveWindow` while `start < chromosome_length − aveWindow` (the last window is not included).
- Value: Pearson correlation coefficient between `numC/(numC+numT)` of the window-1 end and of the window-2 end, over all usable read pairs of that window pair.
- `NaN` is written when fewer than 2 usable pairs exist.

**Behaviour notes**

- Read pairs: records overlapping each window are collected (flag filter, name toggle). A pair is every read name present in both windows. A pair is usable when both ends have at least one CpG call and at least `-minCgNum` CpGs.
- The methylation pattern is fixed to `CG` and the conversion-filter pattern is fixed to `CH`. There are no pattern options.
- Randomness: the control mode uses an unseeded `MersenneTwister`; there is no seed option.
- Run time: every window pair issues two BAM region queries of `-aveWindow` bp and holds all their records in memory. Nothing is cached between window pairs, so each window is read `2 × windows` times. The output is flushed after each row.

**Pitfalls**

- GCH (or any non-CpG context) cannot be requested. For a NOMe-HiC BAM the `WCH` filter cannot be selected either; the only filter-related switch is `-turnOffBisulfiteFilter`.
- Diagonal entries (window against itself) are built from records paired with themselves: a name with both mates in the window is removed by the toggle, every remaining record is found in both collections, and there is no same-alignment check.
- The same self-pairing happens for a record that overlaps the shared boundary of two neighbouring windows when its mate is in neither window.
- Without `-useBadMate` long-range pairs that are not flagged proper give `NaN` (unit test, `global` scenario).
- If `chr` is not in the sizes file, or the chromosome is not longer than `-aveWindow`, the output is empty.

**Minimal example**

```
java -Xmx20G -cp "/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorGlobal CpG_global.chr22 genome.chrom.sizes 22 sample.calmd.bam -useBadMate
```

This measures CpG methylation with the `CH` conversion filter. A GCH/WCH invocation as in README Analysis 4 is not possible with this class.

---

## `MethyCorSepAlleleSC`

**Purpose** — From per-cell BAMs and a VCF of heterozygous SNPs it builds, for every CpG seen on a SNP-overlapping read, a vector with one methylation value (0 or 1) per cell and allele, and then writes two CpG × CpG matrices: Pearson correlations and Poisson-binomial scores between those vectors. Choose it only for small regions or small SNP sets, because the matrices grow with the square of the number of CpGs. It measures `CG` only.

**Command**

```
java -Xmx20G -cp "<repo>/target/bisulfitehic-0.38-jar-with-dependencies.jar:<repo>/lib/ccInference-0.19-slim.jar:<repo>/lib/gatk-package-distribution-3.3.jar" \
  main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorSepAlleleSC [opts] outputPrefix snp.vcf.gz -bams cell1.bam -bams cell2.bam ...
```

Positional arguments:
1. `outputPrefix` (`arguments.get(0)`).
2. `snp_vcf.gz` — VCF of SNPs (`arguments.get(1)`).

`-bams` is required (at least one); otherwise a `CmdLineException` is thrown.

**Inputs**

- VCF: bgzip-compressed with a tabix index next to it (the tool stops with a "please use tabix to index it" error otherwise). The whole file is read. Records used: unfiltered, not indel, not mixed, single-position, and heterozygous with an unfiltered genotype in the **first** sample (index 0 is hard-coded). Allele A is the first genotype allele, allele B the second. Phasing is not checked.
- BAMs (`-bams`): one per cell, coordinate-sorted and indexed, with `MD` tags, duplicate-marked. Chromosome names must match the VCF.

**Options**

| Option | Default | What it does in the code |
|---|---|---|
| `-minBaseQ` | 5 | Minimum base quality for a CpG call and for the SNP base. |
| `-bams` | null | Required, repeatable. One BAM per cell. |
| `-minCellNum` | 2 | Parsed but never used. No effect. |
| `-useBadMate` | false | Keeps records without the proper-pair flag. |
| `-minMapQ` | 20 | Minimum mapping quality. |
| `-turnOffBisulfiteFilter` | false | Skips the conversion filter (fixed to `CH`, position 1). |
| `-skipProcessFirstRow` | false | Parsed but never used. No effect. |
| `-h` | false | No effect of its own (usage is printed only when fewer than 2 tokens are given). |

**Output**

Two gzip files, tab-separated, each with one header line. Let N be the number of distinct CpG positions collected. CpGs appear in `HashMap` order (not sorted by coordinate), in the same order in rows and columns and in both files.

`<outputPrefix>.methyCor.txt.gz`

| Col | Name | Meaning |
|---|---|---|
| 1 | `chr` | Chromosome of the row CpG. |
| 2 | `pos` | Position of the row CpG (see Behaviour notes for how it is computed). |
| 3 … N+2 | `<i>_<pos>` | Header label is the 0-based CpG index, an underscore and the CpG position. Value: Pearson correlation between the row CpG's vector and the column CpG's vector. Written only for column index > row index; the diagonal and the lower triangle are `NaN`. |

`<outputPrefix>.poissonBinomial.txt.gz`

- Same header, same row and column layout, `NaN` on the diagonal and in the lower triangle.
- Value for row CpG r and column CpG c: `get_cdf(k)` of a Poisson-binomial distribution whose success probabilities are the entries of r's vector, where k is the number of methylated entries in c's vector (`NaN` entries count as 0).
- Reading `PoissonBinomial2`: when r's vector has no `NaN`, the value is 1.0 if k is at least the number of methylated entries of r and 0.0 otherwise. When r's vector contains `NaN`, the value is 0.0 if k is smaller than the number of methylated entries of r and `NaN` otherwise. (Derived from the source, not from a run.)

**Behaviour notes**

- Vector layout: length `2 × number of BAMs`; cell i (0-based order of `-bams`) uses slot `2i` for allele A and slot `2i+1` for allele B. Slots start as `NaN`.
- For every BAM and every SNP, records overlapping the SNP are read. Unpaired records are always dropped. For C/T SNPs only G-to-A reads are used; for A/G SNPs only C-to-T reads are used. The SNP base must have quality at least `-minBaseQ`. The read is assigned to allele A or B in bisulfite space.
- Every CpG call on that read (value 1 = methylated, 0 = unmethylated) is written into the read's cell/allele slot of that CpG. A later read overwrites the value of an earlier read for the same CpG, cell and allele; values are not averaged.
- CpG position = alignment start of the read + offset of the CpG's cytosine (forward strand) within the soft-clipped read. Both read strands report the cytosine position of the CpG. The offset is an index into the clipped read, so the code does not correct for insertions or deletions in the CIGAR.
- A CpG whose offset equals the SNP offset is skipped.
- Vectors are passed to the correlation with their `NaN` entries.
- Not verified: that the Pearson value is `NaN` whenever either vector has a `NaN` slot, i.e. unless both CpGs have a call for both alleles in every cell (expected from commons-math, not run).
- No randomness.
- Memory and run time: all SNPs of the VCF and all CpG vectors are held in memory; each BAM is queried once per SNP; the output step is N² correlations and N² `PoissonBinomial2` objects.

**Pitfalls**

- CpG only, `CH` conversion filter only; GCH cannot be requested.
- `-minCellNum` and `-skipProcessFirstRow` do nothing.
- There is no `-sample` option; the first VCF sample is always used.
- Only CpGs on the same read as the SNP base are collected; the mate is not used.
- The final log line names the wrong class (`MethyCorAcrossHiccupsSC's running time`).
- Not verified: behaviour when the VCF header has no contig lines (the sequence dictionary is taken from the VCF header).

**Minimal example**

```
java -Xmx20G -cp "/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/path/to/bisulfitehic/lib/ccInference-0.19-slim.jar:/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorSepAlleleSC sepAllele.region phased_snps.region.vcf.gz -bams cell1.calmd.bam -bams cell2.calmd.bam -useBadMate
```

This measures CpG methylation; the class has no GCH/WCH options. No README or test invocation exists for this class.

---

## Known limitations

Found while reading the source for this guide (November 2026). They are described here, not fixed.

All five tools

1. `-h` does nothing; usage is printed only when too few tokens are given. A missing list option (`-bams`, `-bws`, `-bwCounts`) ends in a stack trace.
2. The supplementary flag (0x800) is never tested, and the read-name toggle adds a name back on a third record.
3. BEDPE and window coordinates are passed to the BAM and bigWig queries as they are, with no conversion from 0-based to 1-based.
4. `-skipProcessFirstRow` skips the first physical line only, also when that line is a `#` comment.
5. `-useUnpaired` without `-useBadMate` stops with `IllegalStateException` on the first unpaired record (`MethyCorAcrossHiccups`, `MethyCorGlobal`).
6. The conversion filter counts pattern matches on the read and on the reference separately, not at the same positions.

`MethyCorAcrossHiccups`

7. `-randomSampleReads` never draws a random subset. With `-randomHiccupsInput` it is the number of shuffle rounds, so each pair is written up to that many times; with `-onlyShort` it stops after N+1 rows per BEDPE row.
8. Output columns 2-3 and 5-6 are the anchor coordinates from the BEDPE row, not read positions.
9. In the `-vcf` output the difference columns are allele B minus allele A, times 100. Phasing is not required: the code tests only that the genotype is heterozygous.
10. Options ignored in some modes: `-outputCount` (random, short, ASM), `-removeReadnameSuffix` and `-minReads` (ASM), `-randomSampleReads` (default mode, `-randomNotLinked`), second chromosome of `-diffChr` (short mode).
11. Differences between modes: `-minReads` is "at least" in default and short mode and "more than" in random mode; the `id` column starts at 1 in default mode and 0 elsewhere; `region_name` has 6 fields in default and random mode and 5 in short and ASM mode.
12. In default and random mode a read lying in both anchors can be paired with itself when the anchors overlap.
13. `-onlyShort 0` switches off the anchor distance filter but still runs the long-range mode.
14. Shuffles are unseeded, so random-mode output differs between runs.

`MethyCorGlobal`, `MethyCorEachChrScBw`

15. The last window of the chromosome is always dropped, even when it is complete.
16. `MethyCorGlobal` pairs a record with itself on the matrix diagonal and where a record overlaps a window boundary.
17. `MethyCorGlobal` is fixed to CpG with a `CH` conversion filter; GCH cannot be requested, and the `WCH` filter for NOMe data cannot be selected.
18. `-aveWindow` is always used; the usage text says the default uses no window.
19. `MethyCorEachChrScBw`: `-random` and `-seed` are ignored with `-shuffleWithinCell`.

`MethyCorAcrossHiccupsSC`

20. `-onlyShort` cannot report pairs inside one anchor (the name toggle removes them).
21. Every BAM is reopened for every BEDPE row, which is slow with many cells.

`MethyCorSepAlleleSC`

22. `-minCellNum` and `-skipProcessFirstRow` are declared and never used.
23. CpG coordinates are alignment start plus offset in the clipped read, with no correction for insertions or deletions.
24. A later read overwrites an earlier one at the same CpG; empty slots are passed to the correlation as `NaN`.
25. The test that skips the SNP position uses the cytosine offset for G-to-A reads, so a SNP at the G of a CpG is not skipped there.
26. CpG only, `CH` conversion filter only.
