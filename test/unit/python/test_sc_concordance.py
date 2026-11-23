"""Unit tests for the scripts in src/python/methy_sc_concordance."""
import gzip
import os
import tempfile
import unittest

import numpy as np
from scipy import stats

import helpers


class PoiBinMixin(object):
    """The same checks run on poibin.py and on the copy embedded in call_methy_concord_sc.py."""
    script = None

    def setUp(self):
        self.PoiBin = helpers.load_script(self.script).PoiBin

    def test_pval_scalar_matches_binomial(self):
        pb = self.PoiBin(np.array([0.01] * 50, dtype='float64'))
        self.assertAlmostEqual(pb.pval(5) / stats.binom.sf(4, 50, 0.01), 1.0, places=6)
        self.assertEqual(pb.pval(0), 1)

    def test_pval_accepts_numpy_int32(self):
        pb = self.PoiBin(np.array([0.2] * 20, dtype='float64'))
        self.assertAlmostEqual(pb.pval(np.int32(5)), stats.binom.sf(4, 20, 0.2), places=10)
        self.assertAlmostEqual(pb.pval(np.int64(5)), stats.binom.sf(4, 20, 0.2), places=10)

    def test_pval_list_input(self):
        pb = self.PoiBin(np.array([0.2] * 20, dtype='float64'))
        pvals = pb.pval([5, 6])
        self.assertAlmostEqual(pvals[0], stats.binom.sf(4, 20, 0.2), places=10)
        self.assertAlmostEqual(pvals[1], stats.binom.sf(5, 20, 0.2), places=10)

    def test_pval_is_clipped_to_unit_interval(self):
        # float32 probabilities make 1 - cdf(k-1) come out slightly negative in the far tail
        pb = self.PoiBin(np.array([0.01] * 50, dtype='float32'))
        for k in range(0, 51):
            p = pb.pval(k)
            self.assertGreaterEqual(p, 0.0, 'k=%d' % k)
            self.assertLessEqual(p, 1.0, 'k=%d' % k)
        for p in pb.pval(list(range(0, 51))):
            self.assertGreaterEqual(p, 0.0)
            self.assertLessEqual(p, 1.0)


class PoiBinModule(PoiBinMixin, unittest.TestCase):
    script = os.path.join('methy_sc_concordance', 'poibin.py')


class PoiBinEmbedded(PoiBinMixin, unittest.TestCase):
    script = os.path.join('methy_sc_concordance', 'call_methy_concord_sc.py')


class CallMethyConcordSc(unittest.TestCase):
    def setUp(self):
        self.mod = helpers.load_script(os.path.join('methy_sc_concordance', 'call_methy_concord_sc.py'))
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.addCleanup(os.chdir, os.getcwd())   # the script chdir()s into its input directory

    def test_unwrap_interval_accepts_any_contig_name(self):
        self.assertEqual(self.mod.unwrap_interval('chr19:100-200'), ('chr19', '100', '200'))
        self.assertEqual(self.mod.unwrap_interval('chrX:100-200'), ('chrX', '100', '200'))
        self.assertEqual(self.mod.unwrap_interval('19:100-200'), ('19', '100', '200'))
        self.assertEqual(self.mod.unwrap_interval('chrUn_gl000220:5-9'), ('chrUn_gl000220', '5', '9'))

    def test_pvalues_are_accurate_and_within_unit_interval(self):
        in_dir = os.path.join(self.tmp.name, 'cells')
        os.makedirs(in_dir)
        n_cell = 12
        for i in range(n_cell):
            with open(os.path.join(in_dir, 'cell_%02d.txt' % i), 'w') as f:
                # window, -, coverage, methylated count, -, methylation probability
                f.write('chr1:100-600\t.\t1\t%d\t.\t0.01\n' % (1 if i < 10 else 0))
                f.write('chrX:100-600\t.\t1\t%d\t.\t0.5\n' % (1 if i < 6 else 0))
        out = os.path.join(self.tmp.name, 'out.txt.gz')
        with helpers.quiet():
            self.mod.call_methy_concordance([in_dir, 'cell_*.txt', out])
        prob = {'chr1': 0.01, 'chrX': 0.5}
        count = {'chr1': 10, 'chrX': 6}
        seen = set()
        with gzip.open(out, 'rt') as f:
            for line in f:
                fields = line.rstrip('\n').split('\t')
                self.assertEqual(len(fields), 8)
                w1, w2, p = fields[0], fields[3], float(fields[7])
                seen.add((w1, w2))
                exact = stats.binom.sf(count[w2] - 1, n_cell, prob[w1])
                self.assertGreaterEqual(p, 0.0)
                self.assertLessEqual(p, 1.0)
                self.assertAlmostEqual(p, exact, delta=1e-12, msg='%s vs %s' % (w1, w2))
        self.assertEqual(seen, {('chr1', 'chr1'), ('chr1', 'chrX'), ('chrX', 'chr1'), ('chrX', 'chrX')})


class CalConcordScAtRegions(unittest.TestCase):
    N_CELL = 30

    @classmethod
    def setUpClass(cls):
        import pyBigWig
        cls.mod = helpers.load_script(os.path.join('methy_sc_concordance', 'cal_concord_sc_at_regions.py'))
        cls.tmp = tempfile.TemporaryDirectory()
        cls.files = []
        for i in range(cls.N_CELL):
            path = os.path.join(cls.tmp.name, 'cell_%02d.bw' % i)
            bw = pyBigWig.open(path, 'w')
            bw.addHeader([('chr1', 10000)])
            v1 = i / float(cls.N_CELL - 1)
            v2 = v1 + 0.01 * ((i * 7) % 5)      # strongly, not perfectly, correlated
            bw.addEntries(['chr1', 'chr1'], [100, 1000], ends=[200, 1100], values=[v1, v2])
            bw.close()
            cls.files.append(path)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    LINE = 'chr1\t100\t200\tchr1\t1000\t1100\n'

    def test_parametric_pvalue_keeps_precision(self):
        fields = self.mod.process_line(10, 0, self.files, self.LINE).rstrip('\n').split('\t')
        self.assertEqual(len(fields), 12)
        self.assertGreater(float(fields[7]), 0.99)
        self.assertGreater(float(fields[8]), 0.0)     # Pearson p-value, about 1e-40
        self.assertLess(float(fields[8]), 1e-20)
        self.assertEqual(fields[11], str(self.N_CELL))

    def test_permutation_pvalue_is_never_zero(self):
        np.random.seed(1)
        n_perm = 20
        fields = self.mod.process_line(10, n_perm, self.files, self.LINE).rstrip('\n').split('\t')
        self.assertEqual(len(fields), 14)
        # the observed correlation beats every permutation: p = (0 + 1) / (n_perm + 1)
        self.assertAlmostEqual(float(fields[12]), 1.0 / (n_perm + 1), places=5)
        self.assertAlmostEqual(float(fields[13]), 1.0 / (n_perm + 1), places=5)

    def test_unparsable_line_is_not_returned_as_output(self):
        self.assertIsNone(self.mod.process_line(10, 0, self.files, 'chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\n'))


class TestIndepPoibin(unittest.TestCase):
    def test_every_row_has_one_column_per_cpg(self):
        mod = helpers.load_script(os.path.join('methy_sc_concordance', 'test_indep_poibin.py'))
        with tempfile.TemporaryDirectory() as tmp:
            in_gz = os.path.join(tmp, 'in.txt.gz')
            out_gz = os.path.join(tmp, 'out.txt.gz')
            with gzip.open(in_gz, 'wt') as f:
                f.write('chr1\t99\t100\t.\t.\t+\t0\t1\t0\t1\tNA\t1\n')
                f.write('chr1\t199\t200\t.\t.\t+\t1\t1\t0\t0\t1\t0\n')
                f.write('chr1\t4999999\t5000000\t.\t.\t+\t0\t1\t1\t0\t0\t1\n')
            mod.test_indep_poibin([in_gz, out_gz])
            with gzip.open(out_gz, 'rt') as f:
                rows = [line.rstrip('\n').split('\t') for line in f]
        self.assertEqual([len(r) for r in rows], [3, 3, 3])
        # pairs further apart than 3 Mb are NA, everything else is a number
        self.assertEqual(rows[0][2], 'NA')
        self.assertEqual(rows[1][2], 'NA')
        self.assertEqual(rows[2][0], 'NA')
        self.assertEqual(rows[2][1], 'NA')
        for i, j in [(0, 0), (0, 1), (1, 0), (1, 1), (2, 2)]:
            float(rows[i][j])


if __name__ == '__main__':
    unittest.main()
