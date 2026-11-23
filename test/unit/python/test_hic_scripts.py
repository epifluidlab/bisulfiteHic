"""Unit tests for mh_reads_summary.py, interact_density_around_regions.py and sam2juicer_new.py."""
import os
import subprocess
import sys
import tempfile
import unittest

import helpers

PAIRED = 0x1
REVERSE = 0x10
MATE_REVERSE = 0x20
READ1 = 0x40
READ2 = 0x80


class MhReadsSummary(unittest.TestCase):
    def setUp(self):
        self.mod = helpers.load_script('mh_reads_summary.py')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_script(self, records):
        bam = helpers.write_bam(os.path.join(self.tmp.name, 'in.bam'), [('chr1', 1000000), ('chr2', 1000000)], records)
        out = os.path.join(self.tmp.name, 'summary.txt')
        self.mod.mh_reads_summary(bam, out)
        summary = {}
        with open(out) as f:
            for line in f:
                fields = line.rstrip('\n').split('\t')
                summary[fields[0]] = fields[1:]
        return summary

    def test_long_range_cis_with_zero_tlen_is_counted(self):
        records = [
            # same-strand pair 50 kb apart: Bhmem writes TLEN=0 for it
            dict(name='a', flag=PAIRED | READ1, ref='chr1', pos=1000, mate_ref='chr1', mate_pos=51000, tlen=0),
            dict(name='a', flag=PAIRED | READ2, ref='chr1', pos=51000, mate_ref='chr1', mate_pos=1000, tlen=0),
            # ordinary forward-reverse pair, 300 bp apart
            dict(name='b', flag=PAIRED | READ1 | MATE_REVERSE, ref='chr1', pos=2000, mate_ref='chr1', mate_pos=2300, tlen=310),
            dict(name='b', flag=PAIRED | READ2 | REVERSE, ref='chr1', pos=2300, mate_ref='chr1', mate_pos=2000, tlen=-310),
            # trans pair
            dict(name='c', flag=PAIRED | READ1, ref='chr1', pos=3000, mate_ref='chr2', mate_pos=3000, tlen=0),
            dict(name='c', flag=PAIRED | READ2, ref='chr2', pos=3000, mate_ref='chr1', mate_pos=3000, tlen=0),
        ]
        s = self.run_script(records)
        for q in ('1', '30'):
            self.assertEqual(s['UniqMappedNoPcrMinMapQ%sTrans:' % q][0], '2')
            self.assertEqual(s['UniqMappedNoPcrMinMapQ%sCis:' % q][0], '4')
            self.assertEqual(s['UniqMappedNoPcrMinMapQ%sCisMore1kb:' % q][0], '2')
            self.assertEqual(s['UniqMappedNoPcrMinMapQ%sCisMore20kb:' % q][0], '2')
        for key, fields in s.items():
            self.assertTrue(fields[1].endswith('%'), key)

    def test_empty_input_does_not_divide_by_zero(self):
        s = self.run_script([])
        self.assertEqual(s['TotalReads:'][0], '0')
        self.assertEqual(s['UniqMappedNoPcrMinMapQ30CisMore20kb:'][0], '0')


class InteractDensity(unittest.TestCase):
    REFS = [('chr1', 20000000), ('chr2', 1000000)]

    def setUp(self):
        self.mod = helpers.load_script('interact_density_around_regions.py')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_script(self, records, from_region, plot_center, plot_bin_size, plot_boundary, min_insert_size=20000):
        bam = helpers.write_bam(os.path.join(self.tmp.name, 'in.bam'), self.REFS, records, sort_index=True)
        out = os.path.join(self.tmp.name, 'density.txt')
        with helpers.quiet():
            self.mod.interact_density_around_regions(bam, out, from_region, plot_center, plot_bin_size,
                                                     plot_boundary, 30, min_insert_size)
        with open(out) as f:
            labels = [int(x) for x in f.readline().rstrip('\n').rstrip('\t').split('\t')]
            counts = [int(x) for x in f.readline().rstrip('\n').rstrip('\t').split('\t')]
        return labels, counts

    def trans_reads(self, center, offsets):
        return [dict(name='r%d' % i, flag=PAIRED | READ1, ref='chr2', pos=2000 + i, mate_ref='chr1',
                     mate_pos=center + off) for i, off in enumerate(offsets)]

    def test_bins_left_of_centre_use_floor(self):
        center = 10000000
        offsets = [-1, 1, -50000, -49999, -5000000, -4999999, 4999999]
        labels, counts = self.run_script(self.trans_reads(center, offsets), 'chr2:1000-5000',
                                         'chr1:%d-%d' % (center, center), 50000, 5000000)
        self.assertEqual(len(counts), 200)
        self.assertEqual(labels[100], 0)
        expected = [0] * 200
        expected[0] = 2     # -5,000,000 and -4,999,999
        expected[99] = 3    # -50,000, -49,999 and -1  (bin [-50 kb, 0))
        expected[100] = 1   # +1                       (bin [0, +50 kb))
        expected[199] = 1   # +4,999,999
        self.assertEqual(counts, expected)

    def test_bin_size_not_dividing_boundary_does_not_crash(self):
        center = 10000000
        offsets = [10, -10, 99, -100]
        labels, counts = self.run_script(self.trans_reads(center, offsets), 'chr2:1000-5000',
                                         'chr1:%d-%d' % (center, center), 30, 100)
        # 6 bins of 30 bp cover [-90, +90); mates outside that range are not counted
        self.assertEqual(counts, [0, 0, 1, 1, 0, 0])

    def test_min_insert_size_uses_mate_positions_not_tlen(self):
        records = [
            # same-strand cis pair only 500 bp apart, TLEN=0 as written by Bhmem: must be filtered
            dict(name='short', flag=PAIRED | READ1, ref='chr1', pos=2600, mate_ref='chr1', mate_pos=3100, tlen=0),
            # same-strand cis pair ~97 kb apart, TLEN=0: must be kept
            dict(name='long', flag=PAIRED | READ1, ref='chr1', pos=100500, mate_ref='chr1', mate_pos=3150, tlen=0),
        ]
        labels, counts = self.run_script(records, 'chr1:2500-101000', 'chr1:3000-3000', 100, 1000)
        self.assertEqual(sum(counts), 1)
        self.assertEqual(counts[labels.index(1)], 1)   # offset +150 -> bin [100, 200)


class Sam2JuicerNew(unittest.TestCase):
    HEADER = ['@HD\tVN:1.6\tSO:queryname', '@SQ\tSN:chr1\tLN:100000', '@SQ\tSN:chr2\tLN:100000',
              '@SQ\tSN:MT\tLN:16569']

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_script(self, sam_lines, frag_lines):
        sam = os.path.join(self.tmp.name, 'in.sam')
        frag = os.path.join(self.tmp.name, 'frag.txt')
        with open(sam, 'w') as f:
            f.write('\n'.join(self.HEADER + sam_lines) + '\n')
        with open(frag, 'w') as f:
            f.write('\n'.join(frag_lines) + '\n')
        cmd = [sys.executable, '-B', os.path.join(helpers.SRC, 'sam2juicer_new.py'), '-s', sam, '-f', frag]
        return subprocess.run(cmd, cwd=self.tmp.name, env=helpers.script_env(self.tmp.name),
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    def test_mt_contig_is_found_in_fragment_file(self):
        sam = ['r1\t65\tMT\t101\t60\t10M\tchr1\t501\t0\tAAAAAAAAAA\tIIIIIIIIII',
               'r1\t129\tchr1\t501\t60\t10M\tMT\t101\t0\tCCCCCCCCCC\tIIIIIIIIII']
        res = self.run_script(sam, ['chr1 1000 5000 100000', 'MT 50 500 16569'])
        self.assertEqual(res.returncode, 0, res.stderr)
        lines = res.stdout.strip().split('\n')
        self.assertEqual(len(lines), 1)
        fields = lines[0].split(' ')
        self.assertEqual(fields[1], 'chr1')
        self.assertEqual(fields[5], 'MT')
        self.assertEqual(fields[7], '1')   # position 100 falls in the 2nd MT fragment (50..500)

    def test_missing_contig_is_reported_by_name(self):
        sam = ['r1\t65\tchr1\t101\t60\t10M\tchr2\t501\t0\tAAAAAAAAAA\tIIIIIIIIII',
               'r1\t129\tchr2\t501\t60\t10M\tchr1\t101\t0\tCCCCCCCCCC\tIIIIIIIIII']
        res = self.run_script(sam, ['chr1 1000 5000 100000'])
        self.assertNotEqual(res.returncode, 0)
        self.assertIn('chr2', res.stderr)
        self.assertIn('missing from the fragment file', res.stderr)
        self.assertNotIn('Traceback', res.stderr)


if __name__ == '__main__':
    unittest.main()
