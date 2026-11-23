"""Shared helpers for the unit tests of the scripts in src/python.

Run all tests from the repository root with:
    python -B -m unittest discover -s test/unit/python
(needs numpy, scipy, pysam, bx-python, pyBigWig, tqdm).
"""
import contextlib
import importlib.util
import io
import os
import re
import sys
import types

sys.dont_write_bytecode = True  # do not litter src/python with __pycache__

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
SRC = os.path.join(REPO, 'src', 'python')

# betabinom_test_mgch.py does "from fastbetabino import *"; the package is optional
# in the test environment, so a tiny stand-in is used when it is not installed.
HAVE_REAL_FASTBETABINO = importlib.util.find_spec('fastbetabino') is not None
_FASTBETABINO_STUB = "def fit_alpha_beta(imps, clicks, *args, **kwargs):\n    return (1.0, 1.0)\n"


def load_script(relpath):
    """Import a script under src/python by path (they are not a package)."""
    if not HAVE_REAL_FASTBETABINO and 'fastbetabino' not in sys.modules:
        stub = types.ModuleType('fastbetabino')
        exec(_FASTBETABINO_STUB, stub.__dict__)
        sys.modules['fastbetabino'] = stub
    path = os.path.join(SRC, relpath)
    name = 'bhic_' + re.sub(r'\W', '_', relpath)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def script_env(tmpdir):
    """Environment for running a script as a subprocess."""
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    if not HAVE_REAL_FASTBETABINO:
        stub_dir = os.path.join(tmpdir, 'stubs')
        os.makedirs(stub_dir, exist_ok=True)
        with open(os.path.join(stub_dir, 'fastbetabino.py'), 'w') as f:
            f.write(_FASTBETABINO_STUB)
        env['PYTHONPATH'] = stub_dir + os.pathsep + env.get('PYTHONPATH', '')
    return env


def quiet():
    """Silence the progress prints of the scripts."""
    return contextlib.redirect_stdout(io.StringIO())


def bed_row(chrom, start, m, c):
    """One line of the '6plus2' bed: bed6 + methylation (%) + coverage."""
    return '\t'.join([chrom, str(start), str(start + 1), '.', '0', '+', str(m), str(c)])


def write_tabix_bed(path, rows):
    """Write rows (already sorted) to path, bgzip + tabix it, return the .gz path."""
    import pysam
    with open(path, 'w') as f:
        for row in rows:
            f.write(row + '\n')
    pysam.tabix_index(path, preset='bed', force=True)
    return path + '.gz'


def write_bam(path, references, records, sort_index=False):
    """Write a tiny BAM.

    references: list of (name, length)
    records: dicts with name, flag, ref, pos, mate_ref, mate_pos and optional mapq, tlen
    """
    import pysam
    header = {'HD': {'VN': '1.6', 'SO': 'unsorted'},
              'SQ': [{'SN': n, 'LN': l} for n, l in references]}
    ref_ids = {n: i for i, (n, _) in enumerate(references)}
    raw = path + '.unsorted.bam' if sort_index else path
    with pysam.AlignmentFile(raw, 'wb', header=header) as out:
        for rec in records:
            a = pysam.AlignedSegment(out.header)
            a.query_name = rec['name']
            a.flag = rec['flag']
            a.reference_id = ref_ids[rec['ref']]
            a.reference_start = rec['pos']
            a.mapping_quality = rec.get('mapq', 60)
            a.cigarstring = '10M'
            a.query_sequence = 'A' * 10
            a.query_qualities = pysam.qualitystring_to_array('I' * 10)
            a.next_reference_id = ref_ids[rec['mate_ref']]
            a.next_reference_start = rec['mate_pos']
            a.template_length = rec.get('tlen', 0)
            out.write(a)
    if sort_index:
        pysam.sort('-o', path, raw)
        pysam.index(path)
        os.remove(raw)
    return path
