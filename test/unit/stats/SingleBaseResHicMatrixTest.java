import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import main.java.edu.mit.compbio.bisulfitehic.utils.SingleBaseResHicMatrix;
import main.java.edu.mit.compbio.bisulfitehic.utils.SingleBaseResHicMatrixV2;

/**
 * End-to-end test of SingleBaseResHicMatrix and SingleBaseResHicMatrixV2 on a synthetic indexed BAM.
 *
 * usage: SingleBaseResHicMatrixTest v1|v2
 *
 * CpGs (C at the 1-based position; the CpG list uses that position as the BED start, which is what
 * both tools look up) at 1000, 5000, 9000, 13000 and these read pairs:
 *   P x3: 1000 meth   - 5000 meth                       -> 1000/5000  : 3 concordant
 *   Q x4: 1000 - 9000, 2 meth/meth + 2 unmeth/unmeth    -> 1000/9000  : 4 concordant
 *   D x2: 1000 meth (mate with a 2 bp deletion upstream of the CpG) - 9000 unmeth
 *                                                       -> 1000/9000  : 2 discordant
 *   R x3: 5000 meth   - 13000 meth                      -> 5000/13000 : 3 concordant
 * No read pair links 1000/13000, 5000/9000 or 9000/13000, so those pairs must not be written.
 */
public class SingleBaseResHicMatrixTest {
	static final int[] CPGS = { 1000, 5000, 9000, 13000 };

	static boolean[] state(int cpg, boolean meth) {
		boolean[] m = new boolean[CPGS.length];
		for (int k = 0; k < CPGS.length; k++) {
			m[k] = CPGS[k] == cpg ? meth : true;
		}
		return m;
	}

	static void pair(List<SAMRecord> recs, SAMFileHeader h, String name, int cpg1, boolean meth1, int cpg2, boolean meth2) {
		recs.add(StatsSynth.mate(h, name, true, cpg1, CPGS, state(cpg1, meth1), 0, 'A', cpg2));
		recs.add(StatsSynth.mate(h, name, false, cpg2, CPGS, state(cpg2, meth2), 0, 'A', cpg1));
	}

	static void run(boolean v2, String[] args) throws Exception {
		if (v2) {
			SingleBaseResHicMatrixV2.main(args);
		} else {
			SingleBaseResHicMatrix.main(args);
		}
	}

	/** data rows keyed by "start1-start2"; checks that every row has as many columns as the header */
	static Map<String, String[]> rows(String what, String prefix) throws Exception {
		Map<String, String[]> rows = new HashMap<String, String[]>();
		int headerCols = -1;
		for (String line : T.readGz(prefix + ".singleBaseResHicMatrix.bed.gz")) {
			String[] tok = line.split("\t", -1);
			if (line.startsWith("#")) {
				headerCols = tok.length;
				continue;
			}
			T.check(tok.length == headerCols, what + ": row has " + tok.length + " columns, header has " + headerCols);
			rows.put(tok[1] + "-" + tok[3], tok);
		}
		return rows;
	}

	static void counts(String what, Map<String, String[]> rows, String key, double concord, double discord) {
		String[] tok = rows.get(key);
		T.check(tok != null, what + ": CpG pair " + key + " is reported");
		if (tok != null) {
			T.check(Double.parseDouble(tok[6]) == concord && Double.parseDouble(tok[7]) == discord,
					what + ": CpG pair " + key + " concordant/discordant (got " + tok[6] + "/" + tok[7] + ", want " + concord + "/" + discord + ")");
		}
	}

	static void keys(String what, Map<String, String[]> rows, String... want) {
		TreeSet<String> w = new TreeSet<String>();
		for (String k : want) {
			w.add(k);
		}
		T.check(new TreeSet<String>(rows.keySet()).equals(w), what + ": reported CpG pairs (got " + new TreeSet<String>(rows.keySet()) + ", want " + w + ")");
	}

	public static void main(String[] args) throws Exception {
		boolean v2 = args.length > 0 && args[0].equals("v2");
		String tag = v2 ? "V2" : "V1";
		File dir = T.tmpDir("SingleBaseResHicMatrixTest_" + tag);

		SAMFileHeader h = StatsSynth.header(SAMFileHeader.SortOrder.coordinate);
		List<SAMRecord> recs = new ArrayList<SAMRecord>();
		for (int i = 0; i < 3; i++) {
			pair(recs, h, "P" + i, 1000, true, 5000, true);
			pair(recs, h, "R" + i, 5000, true, 13000, true);
		}
		for (int i = 0; i < 4; i++) {
			pair(recs, h, "Q" + i, 1000, i < 2, 9000, i < 2);
		}
		for (int i = 0; i < 2; i++) {
			// mate 1: 988-992 matched, 993-994 deleted, 995-1029 matched: the CpG at 1000 is read offset 10, reference offset 12
			recs.add(StatsSynth.mate(h, "D" + i, true, 988, 5, 2, 35, CPGS, state(1000, true), 0, 'A', 8990));
			recs.add(StatsSynth.mate(h, "D" + i, false, 9000, CPGS, state(9000, false), 0, 'A', 998));
		}
		File bam = new File(dir, "reads.bam");
		StatsSynth.writeBam(bam, h, recs);

		File cpg = new File(dir, "cpg.bed");
		StringBuilder bed = new StringBuilder();
		for (int c : CPGS) {
			bed.append("chr1\t").append(c).append("\t").append(c + 1).append("\n");
		}
		T.write(cpg, bed.toString());

		// 1. default options
		String what = tag + " default";
		try {
			String prefix = new File(dir, "default").getPath();
			run(v2, new String[] { "-turnOffBisulfiteFilter", prefix, cpg.getPath(), bam.getPath() });
			Map<String, String[]> rows = rows(what, prefix);
			keys(what, rows, "1000-5000", "1000-9000", "5000-13000");
			counts(what, rows, "1000-5000", 3, 0);
			counts(what + " (read with a deletion upstream of the CpG)", rows, "1000-9000", 4, 2);
			counts(what + " (last CpG of the chromosome)", rows, "5000-13000", 3, 0);
		} catch (Throwable e) {
			e.printStackTrace();
			T.check(false, what + " threw " + e);
		}

		// 2. -permutations 0: no expectation columns, but no crash
		what = tag + " -permutations 0";
		try {
			String prefix = new File(dir, "noperm").getPath();
			run(v2, new String[] { "-turnOffBisulfiteFilter", "-permutations", "0", prefix, cpg.getPath(), bam.getPath() });
			Map<String, String[]> rows = rows(what, prefix);
			counts(what, rows, "1000-5000", 3, 0);
			String[] tok = rows.get("1000-5000");
			if (tok != null) {
				T.check(tok[tok.length - 1].equals("NaN"), what + ": permutation p-value is NaN (got " + tok[tok.length - 1] + ")");
			}
		} catch (Throwable e) {
			e.printStackTrace();
			T.check(false, what + " threw " + e);
		}

		// 3. -bedpe with a header row and -skipProcessFirstRow: only the pair inside the two anchors
		what = tag + " -bedpe -skipProcessFirstRow";
		try {
			File bedpe = new File(dir, "loops.bedpe");
			T.write(bedpe, "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\nchr1\t900\t1100\tchr1\t8900\t9100\n");
			String prefix = new File(dir, "bedpe").getPath();
			run(v2, new String[] { "-turnOffBisulfiteFilter", "-skipProcessFirstRow", "-bedpe", bedpe.getPath(), prefix, cpg.getPath(), bam.getPath() });
			Map<String, String[]> rows = rows(what, prefix);
			keys(what, rows, "1000-9000");
			String[] tok = rows.get("1000-9000");
			if (tok != null) {
				T.check(tok[5].equals("chr1:900:1100:8900:9100"), what + ": name is the bedpe anchor pair (got " + tok[5] + ")");
			}
		} catch (Throwable e) {
			e.printStackTrace();
			T.check(false, what + " threw " + e);
		}

		// 4. -regions (V1 only): a pair is kept when either CpG overlaps a region
		if (!v2) {
			what = tag + " -regions";
			try {
				File regions = new File(dir, "regions.bed");
				T.write(regions, "chr1\t12990\t13010\n");
				String prefix = new File(dir, "regions").getPath();
				run(v2, new String[] { "-turnOffBisulfiteFilter", "-regions", regions.getPath(), prefix, cpg.getPath(), bam.getPath() });
				keys(what, rows(what, prefix), "5000-13000");
			} catch (Throwable e) {
				e.printStackTrace();
				T.check(false, what + " threw " + e);
			}
		}

		T.finish("SingleBaseResHicMatrixTest " + tag);
	}
}
