#!/bin/bash
# Run every unit-test group under test/unit/ against the CURRENT sources.
# Nothing is written into the repository (each group compiles into a temp directory).
#
# Usage:  bash test/unit/run_all.sh
# Env:    JAVA_HOME      JDK 9+ for the Java groups (default: Quest OpenJDK 11)
#         PYSAM_PYTHON   python with numpy/scipy/pysam/bx-python for the Python group
#                        (default: the scnomehic conda env on Quest)
# Exit status 0 only if every group passes.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
PYSAM_PYTHON="${PYSAM_PYTHON:-/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/python}"

failed=""
run_group() {  # run_group <name> <command...>
    local name="$1"
    shift
    echo "################ $name"
    if "$@"; then echo ">>> $name: PASS"; else echo ">>> $name: FAIL"; failed="$failed $name"; fi
    echo
}

run_group concordance bash "$HERE/concordance/run.sh"
run_group allele      bash "$HERE/allele/run.sh"
run_group stats       bash "$HERE/stats/run.sh"
run_group perl        bash "$HERE/perl/run.sh"
run_group python      bash -c "cd '$ROOT' && '$PYSAM_PYTHON' -B -m unittest discover -s test/unit/python"

if [ -z "$failed" ]; then
    echo "ALL UNIT TEST GROUPS PASSED"
    exit 0
fi
echo "FAILED GROUPS:$failed"
exit 1
