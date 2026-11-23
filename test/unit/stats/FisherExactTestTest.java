import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;

import main.java.edu.mit.compbio.bisulfitehic.utils.FisherExactTest;

/**
 * FisherExactTest against exact integer arithmetic.
 * The two-sided p-value must include every table that is exactly as probable as the observed one
 * (R's fisher.test definition); comparing doubles with a bare <= dropped such ties.
 */
public class FisherExactTestTest {

	static BigInteger choose(int n, int k) {
		if (k < 0 || k > n) {
			return BigInteger.ZERO;
		}
		BigInteger r = BigInteger.ONE;
		for (int i = 1; i <= k; i++) {
			r = r.multiply(BigInteger.valueOf(n - k + i)).divide(BigInteger.valueOf(i));
		}
		return r;
	}

	/** exact two-sided p: sum of the probabilities of all tables (same margins) no more probable than the observed one */
	static double exactTwoTailed(int a, int b, int c, int d) {
		int r1 = a + b, c1 = a + c, n = a + b + c + d;
		int lo = Math.max(0, c1 - (n - r1)), hi = Math.min(r1, c1);
		BigInteger base = choose(r1, a).multiply(choose(n - r1, c1 - a));
		BigInteger sum = BigInteger.ZERO;
		for (int x = lo; x <= hi; x++) {
			BigInteger w = choose(r1, x).multiply(choose(n - r1, c1 - x));
			if (w.compareTo(base) <= 0) {
				sum = sum.add(w);
			}
		}
		return new BigDecimal(sum).divide(new BigDecimal(choose(n, c1)), MathContext.DECIMAL64).doubleValue();
	}

	public static void main(String[] args) {
		FisherExactTest t = new FisherExactTest(100);

		// point probability and one-sided tails (tea-tasting table: P(a) = 1,16,36,16,1 over 70)
		T.close(t.getP(3, 1, 1, 3), 16.0 / 70, 1e-12, "getP(3,1,1,3)");
		T.close(t.getRightTailedP(3, 1, 1, 3), 17.0 / 70, 1e-12, "getRightTailedP(3,1,1,3)");
		T.close(t.getLeftTailedP(1, 3, 3, 1), 17.0 / 70, 1e-12, "getLeftTailedP(1,3,3,1)");
		T.close(t.getCumlativeP(3, 1, 1, 3), 17.0 / 70, 1e-12, "getCumlativeP(3,1,1,3)");
		T.close(t.getTwoTailedP(3, 1, 1, 3), 34.0 / 70, 1e-12, "getTwoTailedP(3,1,1,3)");
		T.check(Double.isNaN(new FisherExactTest(10).getTwoTailedP(5, 5, 5, 5)), "n > maxSize returns NaN");

		// tables whose opposite tail holds an equally probable table (values from R fisher.test)
		T.close(t.getTwoTailedP(3, 0, 3, 6), 2.0 / 11, 1e-9, "getTwoTailedP(3,0,3,6) tie in the other tail");
		T.close(t.getTwoTailedP(3, 3, 0, 6), 2.0 / 11, 1e-9, "getTwoTailedP(3,3,0,6) tie in the other tail");
		T.close(t.getTwoTailedP(3, 2, 3, 4), 1.0, 1e-9, "getTwoTailedP(3,2,3,4) most probable table has an equal twin");

		// every 2x2 table with n <= 30 against exact arithmetic
		int tables = 0, wrong = 0;
		String firstWrong = "";
		for (int n = 1; n <= 30; n++) {
			for (int a = 0; a <= n; a++) {
				for (int b = 0; b <= n - a; b++) {
					for (int c = 0; c <= n - a - b; c++) {
						int d = n - a - b - c;
						tables++;
						double got = t.getTwoTailedP(a, b, c, d);
						double want = exactTwoTailed(a, b, c, d);
						if (Double.isNaN(got) || Math.abs(got - want) > 1e-9) {
							wrong++;
							if (firstWrong.isEmpty()) {
								firstWrong = "(" + a + "," + b + "," + c + "," + d + ") got " + got + " want " + want;
							}
						}
					}
				}
			}
		}
		T.check(wrong == 0, "two-sided p differs from exact arithmetic for " + wrong + " of " + tables + " tables with n<=30; first: " + firstWrong);

		T.finish("FisherExactTestTest");
	}
}
