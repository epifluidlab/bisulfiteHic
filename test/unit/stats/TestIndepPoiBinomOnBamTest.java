import java.io.File;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import main.java.edu.mit.compbio.bisulfitehic.utils.TestIndepPoiBinomOnBam;

/**
 * End-to-end test of TestIndepPoiBinomOnBam on two synthetic single-cell SAM files and one het SNP.
 *
 * Each cell has one read pair per allele (SNP A/T at 1005 on mate 1). Mate 1 lies in window 0 and
 * carries two CpGs, mate 2 lies in window 2 and carries two CpGs. Methylated CpGs out of 2:
 *
 *                     window 0   window 2   pooled
 *   cell 0 allele A      0          0        0
 *   cell 0 allele B      1          0        0.25
 *   cell 1 allele A      0          1        0.25
 *   cell 1 allele B      1          1        0.5
 *
 * With the four alleles kept apart, var(pooled) = 1/24 and the mean of the two window variances is
 * 1/12, so the window 0 / window 2 score is 0.5. (If cell 0 allele B and cell 1 allele A share a
 * slot there are three groups 0, 0.25, 0.5 in every vector and the score is 1.)
 *
 * In permutation mode the permuted pooled levels are 0, 0.25 or 0.5, so permuted variances reach
 * the observed one in most permutations and -log10(p) stays far below 1.5. (With an integer
 * division every permuted level is 0 and the score is -log10(1/permutations) = 3.3.)
 */
public class TestIndepPoiBinomOnBamTest {
	static final int[] CPGS = { 1010, 1020, 25000, 25010 };

	static void pair(List<SAMRecord> recs, SAMFileHeader h, String name, char snpBase, int meth1, int meth2) {
		// meth1/meth2: number of methylated CpGs (0 or 1) among the two CpGs of mate 1 / mate 2
		boolean[] m = { meth1 > 0, false, meth2 > 0, false };
		recs.add(StatsSynth.mate(h, name, true, 1000, CPGS, m, 1005, snpBase, 25000));
		recs.add(StatsSynth.mate(h, name, false, 25000, CPGS, m, 0, 'A', 1000));
	}

	static String[] row0(String prefix) throws Exception {
		List<String> lines = T.readGz(prefix + ".matrix.txt.gz");
		T.check(lines.size() == 3, "3 windows (got " + lines.size() + " rows)");
		String[] tok = lines.get(0).split("\t", -1);
		T.check(tok.length == 6, "row = chr, start, end + 3 windows (got " + tok.length + " fields: '" + lines.get(0) + "')");
		return tok;
	}

	public static void main(String[] args) throws Exception {
		File dir = T.tmpDir("TestIndepPoiBinomOnBamTest");
		SAMFileHeader h = StatsSynth.header(SAMFileHeader.SortOrder.unsorted);

		List<SAMRecord> cell0 = new ArrayList<SAMRecord>();
		pair(cell0, h, "c0A", 'A', 0, 0);
		pair(cell0, h, "c0B", 'T', 1, 0);
		List<SAMRecord> cell1 = new ArrayList<SAMRecord>();
		pair(cell1, h, "c1A", 'A', 0, 1);
		pair(cell1, h, "c1B", 'T', 1, 1);
		File sam0 = new File(dir, "cell0.sam");
		File sam1 = new File(dir, "cell1.sam");
		StatsSynth.writeSam(sam0, h, cell0);
		StatsSynth.writeSam(sam1, h, cell1);
		File vcf = new File(dir, "snp.vcf.gz");
		StatsSynth.writeHetSnpVcf(vcf, 1005, "A", "T");

		String[] common = { "-turnOffBisulfiteFilter", "-aveWindow", "10000", "-chr", "chr1", "-bams", sam0.getPath(), "-bams", sam1.getPath() };

		try {
			String prefix = new File(dir, "ratio").getPath();
			List<String> a = new ArrayList<String>(java.util.Arrays.asList(common));
			a.add(vcf.getPath());
			a.add(prefix);
			TestIndepPoiBinomOnBam.main(a.toArray(new String[0]));
			String[] tok = row0(prefix);
			if (tok.length == 6) {
				T.check(tok[0].equals("chr1") && tok[1].equals("0") && tok[2].equals("10000"), "row 0 is chr1 0 10000");
				Double v = T.num(tok[5]);
				T.check(v != null && Math.abs(v - 0.5) < 1e-9, "window 0 / window 2 variance ratio over the 4 alleles of 2 cells (got " + tok[5] + ", want 0.5)");
			}
		} catch (Throwable e) {
			e.printStackTrace();
			T.check(false, "TestIndepPoiBinomOnBam threw " + e);
		}

		try {
			String prefix = new File(dir, "perm").getPath();
			List<String> a = new ArrayList<String>(java.util.Arrays.asList(common));
			a.add("-permutation");
			a.add("2000");
			a.add(vcf.getPath());
			a.add(prefix);
			TestIndepPoiBinomOnBam.main(a.toArray(new String[0]));
			String[] tok = row0(prefix);
			if (tok.length == 6) {
				Double v = T.num(tok[5]);
				T.check(v != null && v >= 0 && v < 1.5, "window 0 / window 2 permutation -log10(p) (got " + tok[5] + ", want < 1.5; 3.3 means every permuted level was 0)");
			}
		} catch (Throwable e) {
			e.printStackTrace();
			T.check(false, "TestIndepPoiBinomOnBam -permutation threw " + e);
		}

		T.finish("TestIndepPoiBinomOnBamTest");
	}
}
