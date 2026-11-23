import java.util.Random;

import main.java.edu.mit.compbio.bisulfitehic.utils.PoissonBinomial2;

/**
 * PoissonBinomial2 PMF/CDF against a direct convolution (pmf[k] = pmf[k]*(1-p) + pmf[k-1]*p).
 * Covers an odd number of trials (the reversed PMF had its two middle entries swapped) and
 * probabilities equal to 1 (the odds p/(1-p) are infinite, which gave NaN).
 */
public class PoissonBinomial2Test {

	static double[] reference(double[] p) {
		double[] f = new double[p.length + 1];
		f[0] = 1;
		for (int i = 0; i < p.length; i++) {
			for (int k = i + 1; k >= 0; k--) {
				f[k] = f[k] * (1 - p[i]) + (k > 0 ? f[k - 1] * p[i] : 0);
			}
		}
		return f;
	}

	static void compare(String name, double[] p) {
		double[] want = reference(p);
		PoissonBinomial2 pb = new PoissonBinomial2(p);
		boolean ok = true;
		String detail = "";
		double cum = 0;
		for (int k = 0; k <= p.length; k++) {
			cum += want[k];
			double pmf = pb.get_pmf(k);
			double cdf = pb.get_cdf(k);
			if (Double.isNaN(pmf) || Math.abs(pmf - want[k]) > 1e-9 || Double.isNaN(cdf) || Math.abs(cdf - cum) > 1e-9) {
				if (ok) {
					detail = " k=" + k + " pmf got " + pmf + " want " + want[k] + ", cdf got " + cdf + " want " + cum;
				}
				ok = false;
			}
		}
		T.check(ok, name + detail);
		T.check(Double.isNaN(pb.get_pmf(p.length + 1)), name + ": more successes than trials is NaN");
	}

	public static void main(String[] args) {
		compare("one trial p=0.3", new double[] { 0.3 });
		compare("three trials (odd)", new double[] { 0.1, 0.5, 0.8 });
		compare("four trials (even)", new double[] { 0.1, 0.5, 0.8, 0.25 });
		compare("five trials (odd)", new double[] { 0.2, 0.2, 0.9, 0.4, 0.6 });
		compare("all zero", new double[] { 0, 0 });
		compare("all one", new double[] { 1, 1 });
		compare("0/1 calls as in MethyCorSepAlleleSC", new double[] { 1, 0, 1, 0, 0, 1 });
		compare("mixed with p=1", new double[] { 1.0, 0.0, 1.0, 0.25 });
		compare("single p=1", new double[] { 1.0 });

		Random rnd = new Random(7);
		for (int n = 1; n <= 12; n++) {
			double[] p = new double[n];
			for (int i = 0; i < n; i++) {
				p[i] = rnd.nextDouble();
			}
			compare("random n=" + n, p);
		}

		T.finish("PoissonBinomial2Test");
	}
}
