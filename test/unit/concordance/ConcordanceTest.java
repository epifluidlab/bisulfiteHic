import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.tuple.Triple;
import org.kohsuke.args4j.CmdLineException;

import main.java.edu.mit.compbio.bisulfitehic.utils.BisulfiteHicUtils;
import main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorAcrossHiccups;
import main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorEachChrScBw;
import main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorGlobal;

/**
 * Unit tests for the methylation-concordance classes (BisulfiteHicUtils, MethyCorAcrossHiccups,
 * MethyCorGlobal, MethyCorEachChrScBw). No real genomic data: every BAM is built in memory from a
 * 31 bp synthetic reference segment and written to the work directory.
 *
 * usage: ConcordanceTest <scenario> <workDir>      (one scenario per JVM, see run.sh)
 * exit code 0 = pass, 1 = fail.
 */
public class ConcordanceTest {

	// CpG at 2, 11, 19, 28; CH cytosines at 6, 15, 24
	static final String REF = "ATCGTTCATTACGTTCTTTCGATTCATACGA";
	static int failures = 0;

	public static void main(String[] args) throws Exception {
		String scenario = args[0];
		File work = new File(args[1]);
		work.mkdirs();
		if (scenario.equals("helpers")) helpers();
		else if (scenario.equals("random")) random(work);
		else if (scenario.equals("count")) count(work);
		else if (scenario.equals("diffchr")) diffchr(work);
		else if (scenario.equals("short")) shortDist(work);
		else if (scenario.equals("selfpair")) selfpair(work);
		else if (scenario.equals("unpaired")) unpaired(work);
		else if (scenario.equals("global")) global(work);
		else if (scenario.equals("oddbws")) oddbws(work);
		else throw new IllegalArgumentException("unknown scenario " + scenario);
		if (failures > 0) {
			System.out.println("FAIL " + scenario + " (" + failures + " check(s) failed)");
			System.exit(1);
		}
		System.out.println("PASS " + scenario);
	}

	static void check(boolean ok, String what) {
		if (!ok) {
			failures++;
			System.out.println("  check failed: " + what);
		}
	}

	// ---------------------------------------------------------------- pure helper tests

	static byte[] quals(int n) {
		byte[] q = new byte[n];
		Arrays.fill(q, (byte) 40);
		return q;
	}

	static Triple<Integer, Integer, Integer> general(String ref, String read, String pat, int cPos, boolean neg, boolean second) {
		return BisulfiteHicUtils.readsMethySummaryWithFivePrimeFilterGeneral(ref.getBytes(), read.getBytes(), quals(ref.length()),
				pat.getBytes(), cPos, neg, second, 5, 0);
	}

	static TreeMap<Integer, Integer> generalLoc(String ref, String read, String pat, int cPos, boolean neg, boolean second) {
		return BisulfiteHicUtils.readsMethyAndLocWithFivePrimeFilterGeneral(ref.getBytes(), read.getBytes(), quals(ref.length()),
				pat.getBytes(), cPos, neg, second, 5, 0);
	}

	static void helpers() {
		// positive controls
		check(general("AACGAA", "AACGAA", "CG", 1, false, false).equals(Triple.of(1, 0, 1)), "CG methylated call on C->T read");
		check(general("AACGAA", "AATGAA", "CG", 1, false, false).equals(Triple.of(0, 1, 1)), "CG unmethylated call on C->T read");
		check(general("TACGAA", "TACGAA", "HCG", 2, false, false).equals(Triple.of(1, 0, 1)), "HCG call");
		check(general("AACGAA", "AACAAA", "CG", 1, true, false).equals(Triple.of(0, 1, 1)), "CG unmethylated call on G->A read");

		// finding 9: '0' (inserted base, as written by htsjdk makeReferenceFromAlignment) and N in the
		// reference must not act as wildcards
		check(general("AA0GAA", "AATGAA", "CG", 1, false, false).equals(Triple.of(0, 0, 0)), "inserted T before G must not be a CpG call");
		check(general("AAC0AA", "AACTAA", "CG", 1, false, false).equals(Triple.of(0, 0, 0)), "C followed by an insertion must not be a CpG call");
		check(general("AA00AA", "AATAAA", "CG", 1, false, false).equals(Triple.of(0, 0, 0)), "2 bp insertion must not be a CpG call");
		check(general("AANGAA", "AACGAA", "CG", 1, false, false).equals(Triple.of(0, 0, 0)), "reference N must not be a CpG call");
		check(general("AAC0AA", "AACTAA", "CG", 1, true, false).equals(Triple.of(0, 0, 0)), "insertion on G->A read must not be a call");
		check(generalLoc("AA0GAA", "AATGAA", "CG", 1, false, false).isEmpty(), "per-location: inserted base must not be a call");

		// finding 9: trailing insertion must not make the reference longer than the clipped read
		check(BisulfiteHicUtils.modifyRefSeqByCigar("ACGT00".getBytes(), "4M2I").length == 4, "4M2I reference length");
		check(BisulfiteHicUtils.modifyRefSeqByCigar("00ACGT0000".getBytes(), "2S4M1I3S").length == 4, "2S4M1I3S reference length");
		check(BisulfiteHicUtils.modifyRefSeqByCigar("00ACGT".getBytes(), "2I4M").length == 4, "2I4M reference length (unchanged behaviour)");
		check(new String(BisulfiteHicUtils.modifyRefSeqByCigar("AC0GT".getBytes(), "2M1I2M")).equals("AC0GT"), "inner insertion kept (unchanged behaviour)");
		check(new String(BisulfiteHicUtils.modifyRefSeqByCigar("ACGT".getBytes(), "2M3D2M")).equals("ACGT"), "deletion (unchanged behaviour)");

		// finding 15: offsets must be reference-forward for G->A reads too
		TreeMap<Integer, Integer> fwd = generalLoc("ACGTTTT", "ACGTTTT", "CG", 1, false, false);
		check(fwd.size() == 1 && Integer.valueOf(1).equals(fwd.get(1)), "C->T read: offset of the C (1), got " + fwd);
		TreeMap<Integer, Integer> rev = generalLoc("ACGTTTT", "ACGTTTT", "CG", 1, true, false);
		check(rev.size() == 1 && Integer.valueOf(1).equals(rev.get(2)), "G->A read: reference-forward offset of the G (2), got " + rev);
		TreeMap<Integer, Integer> rev2 = generalLoc("TTCGATT", "TTCAATT", "HCG", 2, false, true);
		check(rev2.size() == 1 && Integer.valueOf(0).equals(rev2.get(3)), "G->A read, HCG: offset 3 unmethylated, got " + rev2);
	}

	// ---------------------------------------------------------------- synthetic BAM helpers

	static SAMFileHeader header() {
		SAMFileHeader h = new SAMFileHeader();
		h.addSequence(new SAMSequenceRecord("chr1", 5000000));
		h.addSequence(new SAMSequenceRecord("chr2", 5000000));
		h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		return h;
	}

	/** bisulfite-convert REF as seen on a C->T strand read (SEQ on the reference forward strand) */
	static String convert(boolean methylated) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < REF.length(); i++) {
			char c = REF.charAt(i);
			if (c == 'C') {
				boolean cpg = i + 1 < REF.length() && REF.charAt(i + 1) == 'G';
				sb.append(cpg && methylated ? 'C' : 'T');
			} else {
				sb.append(c);
			}
		}
		return sb.toString();
	}

	static String md(String ref, String read) {
		StringBuilder sb = new StringBuilder();
		int run = 0;
		for (int i = 0; i < ref.length(); i++) {
			if (ref.charAt(i) == read.charAt(i)) {
				run++;
			} else {
				sb.append(run).append(ref.charAt(i));
				run = 0;
			}
		}
		sb.append(run);
		return sb.toString();
	}

	/**
	 * A C->T strand read: first-of-pair forward, or second-of-pair reverse. mateChr == null makes it unpaired.
	 */
	static SAMRecord rec(SAMFileHeader h, String name, String chr, int start, boolean second, boolean proper,
			String mateChr, int mateStart, int tlen, boolean methylated) {
		SAMRecord r = new SAMRecord(h);
		String read = convert(methylated);
		char[] q = new char[REF.length()];
		Arrays.fill(q, 'I');
		r.setReadName(name);
		r.setReferenceName(chr);
		r.setAlignmentStart(start);
		r.setMappingQuality(60);
		r.setCigarString(REF.length() + "M");
		r.setReadString(read);
		r.setBaseQualityString(new String(q));
		r.setAttribute("MD", md(REF, read));
		if (mateChr != null) {
			r.setReadPairedFlag(true);
			r.setFirstOfPairFlag(!second);
			r.setSecondOfPairFlag(second);
			r.setProperPairFlag(proper);
			r.setMateUnmappedFlag(false);
			r.setMateReferenceName(mateChr);
			r.setMateAlignmentStart(mateStart);
			r.setReadNegativeStrandFlag(second);
			r.setMateNegativeStrandFlag(!second);
			r.setInferredInsertSize(tlen);
		}
		return r;
	}

	/** both mates of a C->T strand pair */
	static void pair(List<SAMRecord> out, SAMFileHeader h, String name, String chr1, int s1, String chr2, int s2,
			boolean proper, int tlen, boolean meth1, boolean meth2) {
		out.add(rec(h, name, chr1, s1, false, proper, chr2, s2, tlen, meth1));
		out.add(rec(h, name, chr2, s2, true, proper, chr1, s1, -tlen, meth2));
	}

	static String writeBam(File work, String name, SAMFileHeader h, List<SAMRecord> recs) {
		File bam = new File(work, name + ".bam");
		SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(h, false, bam);
		for (SAMRecord r : recs) w.addAlignment(r);
		w.close();
		return bam.getPath();
	}

	static String writeText(File work, String name, String content) throws Exception {
		File f = new File(work, name);
		FileWriter fw = new FileWriter(f);
		fw.write(content);
		fw.close();
		return f.getPath();
	}

	static List<String> readGz(String path) throws Exception {
		List<String> lines = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))));
		String line;
		while ((line = br.readLine()) != null) lines.add(line);
		br.close();
		return lines;
	}

	static List<String> runHiccups(File work, String tag, String bedpe, String bam, String... opts) throws Exception {
		String prefix = new File(work, tag).getPath();
		List<String> a = new ArrayList<String>(Arrays.asList(opts));
		a.add(prefix);
		a.add(bedpe);
		a.add(bam);
		new MethyCorAcrossHiccups().doMain(a.toArray(new String[a.size()]));
		return readGz(prefix + ".methyCor.txt.gz");
	}

	static List<SAMRecord> threeLongRangePairs(SAMFileHeader h) {
		List<SAMRecord> recs = new ArrayList<SAMRecord>();
		pair(recs, h, "p1", "chr1", 1500, "chr1", 50500, true, 49031, true, true);
		pair(recs, h, "p2", "chr1", 1600, "chr1", 50600, true, 49031, false, false);
		pair(recs, h, "p3", "chr1", 1700, "chr1", 50700, true, 49031, true, false);
		return recs;
	}

	// ---------------------------------------------------------------- scenarios

	/** finding 2: -randomHiccupsInput with the default -randomSampleReads must still write rows */
	static void random(File work) throws Exception {
		SAMFileHeader h = header();
		String bam = writeBam(work, "random", h, threeLongRangePairs(h));
		String bedpe = writeText(work, "random.bedpe", "chr1\t1000\t2000\tchr1\t50000\t51000\n");
		List<String> lines = runHiccups(work, "random_out", bedpe, bam, "-useGeneralBedPe", "-randomHiccupsInput");
		check(lines.size() == 3, "expected 3 shuffled rows (one round over 3 shared read pairs), got " + lines.size());
	}

	/** finding 5: -outputCount last column (total read pairs in the BAM) must be the same on every row */
	static void count(File work) throws Exception {
		SAMFileHeader h = header();
		String bam = writeBam(work, "count", h, threeLongRangePairs(h));
		String bedpe = writeText(work, "count.bedpe", "chr1\t1000\t2000\tchr1\t50000\t51000\nchr1\t1000\t2000\tchr1\t50000\t51000\n");
		List<String> lines = runHiccups(work, "count_out", bedpe, bam, "-useGeneralBedPe", "-outputCount");
		check(lines.size() == 6, "expected 6 rows, got " + lines.size());
		for (String l : lines) {
			String[] f = l.split("\t");
			check(f[f.length - 1].equals("3"), "total read pair column should be 3, got " + f[f.length - 1]);
		}
	}

	/** finding 8: with -diffChr the anchor distance filter must not compare coordinates of two chromosomes */
	static void diffchr(File work) throws Exception {
		SAMFileHeader h = header();
		List<SAMRecord> recs = new ArrayList<SAMRecord>();
		pair(recs, h, "x1", "chr1", 1500, "chr2", 5500, false, 0, true, true);
		String bam = writeBam(work, "diffchr", h, recs);
		String bedpe = writeText(work, "diffchr.bedpe", "chr1\t1000\t2000\tchr2\t5000\t6000\n");
		List<String> lines = runHiccups(work, "diffchr_out", bedpe, bam, "-useGeneralBedPe", "-diffChr", "-useBadMate");
		check(lines.size() == 1, "expected 1 inter-chromosomal row, got " + lines.size());
	}

	/** finding 10: -onlyShort / -shortLowerBound must use the mate distance, not TLEN (0 for chimeric pairs) */
	static void shortDist(File work) throws Exception {
		SAMFileHeader h = header();
		List<SAMRecord> recs = new ArrayList<SAMRecord>();
		pair(recs, h, "near", "chr1", 1200, "chr1", 1400, true, 0, true, true);     // 200 bp apart, TLEN 0
		pair(recs, h, "far", "chr1", 1500, "chr1", 100500, true, 0, true, true);    // 99 kb apart, TLEN 0
		String bam = writeBam(work, "short", h, recs);
		String bedpe = writeText(work, "short.bedpe", "chr1\t1000\t2000\tchr1\t100000\t101000\n");
		List<String> lines = runHiccups(work, "short_out1", bedpe, bam, "-useGeneralBedPe", "-onlyShort", "1000");
		check(lines.size() == 1, "-onlyShort 1000: expected only the 200 bp pair, got " + lines.size() + " rows");
		lines = runHiccups(work, "short_out2", bedpe, bam, "-useGeneralBedPe", "-onlyShort", "1000", "-shortLowerBound", "100");
		check(lines.size() == 1, "-onlyShort 1000 -shortLowerBound 100: expected the 200 bp pair, got " + lines.size() + " rows");
		lines = runHiccups(work, "short_out3", bedpe, bam, "-useGeneralBedPe", "-onlyShort", "1000", "-shortLowerBound", "300");
		check(lines.size() == 0, "-onlyShort 1000 -shortLowerBound 300: expected no pair, got " + lines.size() + " rows");
	}

	/** finding 6: in -onlyShort mode a read overlapping both anchors must not be paired with itself */
	static void selfpair(File work) throws Exception {
		SAMFileHeader h = header();
		List<SAMRecord> recs = new ArrayList<SAMRecord>();
		pair(recs, h, "near", "chr1", 1200, "chr1", 1400, true, 231, true, true);   // both mates inside the anchor
		pair(recs, h, "edge", "chr1", 1990, "chr1", 2300, true, 341, true, true);   // mate 2 outside the anchor
		String bam = writeBam(work, "selfpair", h, recs);
		String bedpe = writeText(work, "selfpair.bedpe", "chr1\t1000\t2000\tchr1\t1000\t2000\n");
		List<String> lines = runHiccups(work, "selfpair_out", bedpe, bam, "-useGeneralBedPe", "-onlyShort", "1000");
		check(lines.size() == 1, "expected only the real pair, got " + lines.size() + " rows");
	}

	/** finding 4: -useUnpaired reads (no 0x1 flag) must reach the methylation caller */
	static void unpaired(File work) throws Exception {
		SAMFileHeader h = header();
		List<SAMRecord> recs = new ArrayList<SAMRecord>();
		recs.add(rec(h, "u1/1", "chr1", 1500, false, false, null, 0, 0, true));
		recs.add(rec(h, "u1/2", "chr1", 50500, false, false, null, 0, 0, false));
		String bam = writeBam(work, "unpaired", h, recs);
		String bedpe = writeText(work, "unpaired.bedpe", "chr1\t1000\t2000\tchr1\t50000\t51000\n");
		List<String> lines = runHiccups(work, "unpaired_out", bedpe, bam, "-useGeneralBedPe", "-useUnpaired", "-useBadMate", "-removeReadnameSuffix");
		check(lines.size() == 1, "expected 1 row, got " + lines.size());
		if (lines.size() == 1) {
			String[] f = lines.get(0).split("\t");
			check(f[11].equals("1.0") && f[13].equals("0.0"), "methylation of the two ends should be 1.0 / 0.0, got " + f[11] + " / " + f[13]);
		}
		// without -useUnpaired the unpaired reads stay excluded
		lines = runHiccups(work, "unpaired_out2", bedpe, bam, "-useGeneralBedPe", "-useBadMate", "-removeReadnameSuffix");
		check(lines.size() == 0, "without -useUnpaired expected 0 rows, got " + lines.size());
	}

	/** finding 1: MethyCorGlobal -useBadMate must let non-proper pairs through the methylation caller */
	static void global(File work) throws Exception {
		SAMFileHeader h = header();
		List<SAMRecord> recs = new ArrayList<SAMRecord>();
		pair(recs, h, "g1", "chr1", 100000, "chr1", 1500000, false, 0, true, true);
		pair(recs, h, "g2", "chr1", 100100, "chr1", 1500100, false, 0, false, false);
		pair(recs, h, "g3", "chr1", 100200, "chr1", 1500200, false, 0, true, true);
		pair(recs, h, "g4", "chr1", 100300, "chr1", 1500300, false, 0, false, true);
		String bam = writeBam(work, "global", h, recs);
		String sizes = writeText(work, "global.sizes", "chr1\t3000000\n");
		String prefix = new File(work, "global_out").getPath();
		new MethyCorGlobal().doMain(new String[] { "-useBadMate", prefix, sizes, "chr1", bam });
		List<String> lines = readGz(prefix + ".methyCor.txt.gz");
		check(lines.size() == 2, "expected a 2 x 2 window matrix, got " + lines.size() + " rows");
		if (lines.size() == 2) {
			String cell = lines.get(0).split("\t")[1];
			// x = (1,0,1,0), y = (1,0,1,1): Pearson r = 1/sqrt(3)
			check(!cell.equals("NaN") && Math.abs(Double.parseDouble(cell) - 0.5773502691896258) < 1e-9, "window 0 vs window 1 correlation should be 0.57735, got " + cell);
		}
		// without -useBadMate the non-proper pairs stay excluded
		String prefix2 = new File(work, "global_out2").getPath();
		new MethyCorGlobal().doMain(new String[] { prefix2, sizes, "chr1", bam });
		lines = readGz(prefix2 + ".methyCor.txt.gz");
		check(lines.size() == 2 && lines.get(0).split("\t")[1].equals("NaN"), "without -useBadMate window 0 vs window 1 should be NaN");
	}

	/** finding 13: -shuffleWithinCell with an odd number of -bws must be rejected as a usage error */
	static void oddbws(File work) throws Exception {
		String sizes = writeText(work, "oddbws.sizes", "chr1\t3000000\n");
		String prefix = new File(work, "oddbws_out").getPath();
		try {
			new MethyCorEachChrScBw().doMain(new String[] { "-shuffleWithinCell", "-bws", new File(work, "missing.bw").getPath(), prefix, sizes, "chr1" });
			check(false, "expected a CmdLineException for an odd number of -bws");
		} catch (CmdLineException e) {
			// expected
		} catch (Exception e) {
			check(false, "expected a CmdLineException, got " + e);
		}
	}
}
