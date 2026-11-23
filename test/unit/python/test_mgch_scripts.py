"""Unit tests for betabinom_test_mgch.py, binom_test_mgch.py and smooth_mgch.py."""
import gzip
import os
import signal
import subprocess
import sys
import tempfile
import unittest

from scipy import stats
import numpy as np

import helpers


def read_rows(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        return [line.rstrip('\n').split('\t') for line in f if line.strip()]


class BinomTestMgch(unittest.TestCase):
    def setUp(self):
        self.mod = helpers.load_script('binom_test_mgch.py')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_script(self, rows):
        in_bed = os.path.join(self.tmp.name, 'in.bed')
        out_bed = os.path.join(self.tmp.name, 'out.bed')
        with open(in_bed, 'w') as f:
            for row in rows:
                f.write(row + '\n')
        with helpers.quiet():
            self.mod.binom_test_mgch(in_bed, out_bed, min_gch=2, min_point=10, window=100,
                                     max_window=1000, global_acc=0.285)
        return read_rows(out_bed)

    def test_methylated_count_is_rounded_not_truncated(self):
        # 33.33% of 3 reads is 1 methylated read (0.9999 must not be truncated to 0).
        rows = self.run_script([helpers.bed_row('chr1', 950, 50.0, 10),
                                helpers.bed_row('chr1', 1000, 33.33, 3),
                                helpers.bed_row('chr1', 1050, 50.0, 10)])
        # background = 10 / 20 = 0.5, P(X >= 1 | n=3, p=0.5) = 0.875
        self.assertAlmostEqual(float(rows[1][-1]), 0.875, places=6)

    def test_small_pvalue_keeps_precision(self):
        rows = self.run_script([helpers.bed_row('chr1', 950, 10.0, 100),
                                helpers.bed_row('chr1', 1000, 100.0, 30),
                                helpers.bed_row('chr1', 1050, 10.0, 100)])
        # background = 20 / 200 = 0.1, P(X >= 30 | n=30, p=0.1) = 1e-30
        p = float(rows[1][-1])
        self.assertGreater(p, 0.0)
        self.assertLess(abs(p - 1e-30) / 1e-30, 1e-4)


class BetabinomTestMgch(unittest.TestCase):
    # 20 GCH sites 10 bp apart, so the +/-2.5 kb and +/-5 kb windows hold the same sites.
    SPECIAL = {3: (33.33, 3), 7: (57.14, 7), 11: (28.57, 7), 15: (0.0, 20)}
    M_GLOBAL = 16 * 5 + 1 + 4 + 2 + 0   # methylated counts, rounded
    C_GLOBAL = 16 * 10 + 3 + 7 + 7 + 20

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        rows = []
        self.sites = []
        for i in range(20):
            m, c = self.SPECIAL.get(i, (50.0, 10))
            self.sites.append((m, c))
            rows.append(helpers.bed_row('chr1', 1000 + 10 * i, m, c))
        self.in_bed = helpers.write_tabix_bed(os.path.join(self.tmp.name, 'in.bed'), rows)
        self.chrom_size = os.path.join(self.tmp.name, 'chrom.size')
        with open(self.chrom_size, 'w') as f:
            f.write('chr1\t100000\n')

    def run_function(self):
        mod = helpers.load_script('betabinom_test_mgch.py')
        out_bed = os.path.join(self.tmp.name, 'out.bed')
        with helpers.quiet():
            mod.betabinom_test_mgch(self.in_bed, out_bed, self.chrom_size, min_gch=5, min_point=10,
                                    window1=2500, window2=5000, binomial_only=True, two_sided=False)
        return read_rows(out_bed)

    def run_cli(self, name, *options):
        out_bed = os.path.join(self.tmp.name, name)
        cmd = [sys.executable, '-B', os.path.join(helpers.SRC, 'betabinom_test_mgch.py'),
               '--in_bed', self.in_bed, '--out_bed', out_bed, '--chrom_size', self.chrom_size,
               '--min_gch', '5', '--min_point', '10'] + list(options)
        res = subprocess.run(cmd, cwd=self.tmp.name, env=helpers.script_env(self.tmp.name),
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        self.assertEqual(res.returncode, 0, res.stderr)
        return read_rows(out_bed)

    def test_window2_uses_log10_like_window1(self):
        rows = self.run_function()
        self.assertEqual(len(rows), 20)
        for fields in rows:
            # identical background in both windows -> identical -log10(p)
            self.assertAlmostEqual(float(fields[11]), float(fields[12]), places=3)
        # the 0/20 site: P(X <= 0 | n=20, p=87/177)
        expected = -np.log10(stats.binom.cdf(0, 20, (self.M_GLOBAL - 0) / (self.C_GLOBAL - 20)))
        self.assertAlmostEqual(float(rows[15][11]), expected, places=3)
        self.assertAlmostEqual(float(rows[15][12]), expected, places=3)

    def test_methylated_count_is_rounded_not_truncated(self):
        rows = self.run_function()
        for (m, c), fields in zip(self.sites, rows):
            self.assertAlmostEqual(float(fields[10]), m / 100.0 - self.M_GLOBAL / self.C_GLOBAL, places=3)

    def test_boolean_options_parse_false(self):
        less = self.run_cli('less.bed', '--binomial_only', 'True', '--two_sided', 'False')
        two = self.run_cli('two.bed', '--binomial_only', 'True', '--two_sided', 'True')
        self.assertEqual(len(less[15]), 18)
        expected_less = -np.log10(stats.binom.cdf(0, 20, self.M_GLOBAL / (self.C_GLOBAL - 20)))
        self.assertAlmostEqual(float(less[15][11]), expected_less, places=3)
        # "--two_sided False" must not behave like "--two_sided True"
        self.assertNotAlmostEqual(float(less[15][11]), float(two[15][11]), places=2)
        # "--binomial_only False" adds the beta-binomial column
        full = self.run_cli('full.bed', '--binomial_only', 'False', '--two_sided', 'False')
        self.assertEqual(len(full[15]), 19)


class SmoothMgch(unittest.TestCase):
    def setUp(self):
        self.mod = helpers.load_script('smooth_mgch.py')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_script(self, sites, size, s, min_gch, timeout=60):
        rows = [helpers.bed_row('chr1', pos, m, 10) for pos, m in sites]
        in_bed = helpers.write_tabix_bed(os.path.join(self.tmp.name, 'in.bed'), rows)
        chrom_size = os.path.join(self.tmp.name, 'chrom.size')
        with open(chrom_size, 'w') as f:
            f.write('chr1\t%d\n' % size)
        out_bed = os.path.join(self.tmp.name, 'out.bed')

        def on_alarm(signum, frame):
            raise TimeoutError('smooth_mgch did not terminate within %d s' % timeout)

        old = signal.signal(signal.SIGALRM, on_alarm)
        signal.alarm(timeout)
        try:
            with helpers.quiet():
                self.mod.smooth_mgch(in_bed, out_bed, chrom_size, window=1000, max_window=10000,
                                     min_gch=min_gch, column=7, s=s)
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)
        return {int(f[1]): float(f[3]) for f in read_rows(out_bed)}

    def test_smoothing_factor_option_is_used(self):
        # with --s 0 the degree-1 spline interpolates: every site keeps its own value
        base = 1000000
        sites = [(base + off, m) for off, m in
                 [(0, 0.0), (10, 100.0), (25, 0.0), (40, 100.0), (60, 0.0), (100, 100.0), (130, 0.0), (170, 100.0)]]
        out = self.run_script(sites, size=2000000, s=0, min_gch=3)
        for pos, m in sites:
            self.assertAlmostEqual(out[pos], m, places=2, msg='site at %d' % pos)
        # halfway between 0% and 100%
        self.assertAlmostEqual(out[base + 5], 50.0, places=2)

    def test_terminates_on_short_contig_with_few_sites(self):
        sites = [(1000 + 100 * i, 50.0) for i in range(5)]
        out = self.run_script(sites, size=5000, s=1, min_gch=100, timeout=30)
        self.assertIn(1000, out)
        self.assertIn(1400, out)


if __name__ == '__main__':
    unittest.main()
