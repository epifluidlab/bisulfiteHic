import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Minimal assertion helpers shared by the unit tests in this folder (no test framework needed).
 * Every test class has a main() that ends with T.finish(name), which exits non-zero if any check failed.
 */
public class T {
	static int passed = 0;
	static int failed = 0;

	static void check(boolean ok, String what) {
		if (ok) {
			passed++;
		} else {
			failed++;
			System.err.println("FAIL: " + what);
		}
	}

	static void close(double got, double want, double tol, String what) {
		check(!Double.isNaN(got) && Math.abs(got - want) <= tol, what + " (got " + got + ", want " + want + ")");
	}

	static void finish(String name) {
		System.err.println(name + ": " + passed + " checks passed, " + failed + " failed");
		System.exit(failed > 0 ? 1 : 0);
	}

	/** scratch directory for one test: -Dtest.tmp=<dir> (set by run.sh) or java.io.tmpdir */
	static File tmpDir(String name) {
		File d = new File(System.getProperty("test.tmp", System.getProperty("java.io.tmpdir")), name);
		d.mkdirs();
		return d;
	}

	static void write(File f, String text) throws Exception {
		FileWriter w = new FileWriter(f);
		w.write(text);
		w.close();
	}

	static List<String> readGz(String path) throws Exception {
		List<String> lines = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path)), "UTF-8"));
		String line;
		while ((line = br.readLine()) != null) {
			lines.add(line);
		}
		br.close();
		return lines;
	}

	/** "NA"/"NaN" -> NaN, otherwise the number; an unparsable token (e.g. "0.5NA") -> null */
	static Double num(String token) {
		if (token.equals("NA") || token.equals("NaN")) {
			return Double.NaN;
		}
		try {
			return Double.parseDouble(token);
		} catch (NumberFormatException e) {
			return null;
		}
	}
}
