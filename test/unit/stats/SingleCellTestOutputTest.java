import java.io.File;
import java.util.List;

import main.java.edu.mit.compbio.bisulfitehic.utils.SingleCellTest;

/**
 * End-to-end check of the SingleCellTest output matrix on a tiny CpG x cell matrix:
 * one row per CpG, first field is the CpG coordinate as a number, then one tab-separated
 * score (or NA) per CpG. Also checks that -regions restricts the CpGs.
 */
public class SingleCellTestOutputTest {

	public static void main(String[] args) throws Exception {
		File dir = T.tmpDir("SingleCellTestOutputTest");
		File in = new File(dir, "matrix.txt");
		// chr start end + 3 unused columns + 5 cells. 65546 = 0x1000A: written as a char it would be a newline
		T.write(in, "#chr\tstart\tend\ta\tb\tc\tcell1\tcell2\tcell3\tcell4\tcell5\n"
				+ "chr1\t99\t100\t.\t.\t.\t1\t0\t1\t0\t1\n"
				+ "chr1\t199\t200\t.\t.\t.\t1\t0\t1\t0\tNA\n"
				+ "chr1\t299\t300\t.\t.\t.\tNA\tNA\tNA\t0\t1\n"
				+ "chr1\t65545\t65546\t.\t.\t.\t0\t1\t0\t1\t0\n");

		// Pearson correlation (default): <=2 shared cells -> NA
		String[] labels = { "100", "200", "300", "65546" };
		double nan = Double.NaN;
		double[][] want = { { 1, 1, nan, -1 }, { 1, 1, nan, -1 }, { nan, nan, nan, nan }, { -1, -1, nan, 1 } };
		String prefix = new File(dir, "all").getPath();
		List<String> lines = null;
		try {
			SingleCellTest.main(new String[] { in.getPath(), prefix });
			lines = T.readGz(prefix + ".matrix.txt.gz");
		} catch (Throwable e) {
			e.printStackTrace();
			T.check(false, "SingleCellTest threw " + e);
		}
		if (lines != null) {
			T.check(lines.size() == 4, "one output row per CpG (got " + lines.size() + " rows, want 4)");
			for (int i = 0; i < Math.min(4, lines.size()); i++) {
				String[] tok = lines.get(i).split("\t", -1);
				T.check(tok.length == 5, "row " + i + " has label + 4 columns (got " + tok.length + " fields: '" + lines.get(i) + "')");
				T.check(tok[0].equals(labels[i]), "row " + i + " label is the CpG coordinate (got '" + tok[0] + "', want " + labels[i] + ")");
				for (int j = 0; j < 4 && j + 1 < tok.length; j++) {
					Double v = T.num(tok[j + 1]);
					boolean ok = v != null && (Double.isNaN(want[i][j]) ? Double.isNaN(v) : Math.abs(v - want[i][j]) < 1e-9);
					T.check(ok, "row " + i + " column " + j + " (got '" + tok[j + 1] + "', want " + want[i][j] + ")");
				}
			}
		}

		// -regions: only the CpG at 199-200 overlaps chr1:150-250
		File regions = new File(dir, "regions.bed");
		T.write(regions, "chr1\t150\t250\n");
		String prefix2 = new File(dir, "regions").getPath();
		lines = null;
		try {
			SingleCellTest.main(new String[] { "-regions", regions.getPath(), in.getPath(), prefix2 });
			lines = T.readGz(prefix2 + ".matrix.txt.gz");
		} catch (Throwable e) {
			e.printStackTrace();
			T.check(false, "SingleCellTest -regions threw " + e);
		}
		if (lines != null) {
			T.check(lines.size() == 1, "-regions keeps only the overlapping CpG (got " + lines.size() + " rows, want 1)");
			if (lines.size() >= 1) {
				T.check(lines.get(0).split("\t", -1)[0].equals("200"), "-regions keeps the CpG at 200 (row: '" + lines.get(0) + "')");
			}
		}

		T.finish("SingleCellTestOutputTest");
	}
}
