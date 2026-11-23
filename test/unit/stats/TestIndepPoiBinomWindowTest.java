import java.io.File;
import java.util.List;

import main.java.edu.mit.compbio.bisulfitehic.utils.TestIndepPoiBinom;

/**
 * End-to-end check of the TestIndepPoiBinom output shape on a tiny CpG x cell matrix:
 * the result is a square matrix with one row and one column per CpG (aveWindow=1) or per
 * window that holds data (aveWindow>1), including the last window of the file.
 */
public class TestIndepPoiBinomWindowTest {

	static void shape(String what, List<String> lines, int want) {
		T.check(lines.size() == want, what + ": " + want + " rows (got " + lines.size() + ")");
		for (int i = 0; i < lines.size(); i++) {
			// every value is followed by a tab, so drop trailing empty fields
			int n = lines.get(i).split("\t").length;
			T.check(n == want, what + ": row " + i + " has " + want + " columns (got " + n + ")");
		}
	}

	public static void main(String[] args) throws Exception {
		File dir = T.tmpDir("TestIndepPoiBinomWindowTest");
		File in = new File(dir, "matrix.txt");
		// 6 CpGs, two in each of the windows (0,10], (10,20], (20,30]; 4 cells, no missing values
		T.write(in, "chr1\t0\t1\t.\t.\t.\t1\t0\t1\t0\n"
				+ "chr1\t1\t2\t.\t.\t.\t1\t0\t0\t0\n"
				+ "chr1\t10\t11\t.\t.\t.\t1\t1\t1\t0\n"
				+ "chr1\t11\t12\t.\t.\t.\t0\t0\t1\t0\n"
				+ "chr1\t20\t21\t.\t.\t.\t0\t1\t1\t0\n"
				+ "chr1\t21\t22\t.\t.\t.\t1\t1\t0\t0\n");

		String p1 = new File(dir, "single").getPath();
		try {
			TestIndepPoiBinom.main(new String[] { in.getPath(), p1 });
			shape("aveWindow=1", T.readGz(p1 + ".matrix.txt.gz"), 6);
		} catch (Throwable e) {
			e.printStackTrace();
			T.check(false, "TestIndepPoiBinom threw " + e);
		}

		String p2 = new File(dir, "window").getPath();
		try {
			TestIndepPoiBinom.main(new String[] { "-aveWindow", "10", in.getPath(), p2 });
			shape("aveWindow=10 (3 windows with data, last one included)", T.readGz(p2 + ".matrix.txt.gz"), 3);
		} catch (Throwable e) {
			e.printStackTrace();
			T.check(false, "TestIndepPoiBinom -aveWindow 10 threw " + e);
		}

		T.finish("TestIndepPoiBinomWindowTest");
	}
}
