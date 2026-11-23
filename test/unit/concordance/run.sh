#!/bin/bash
# Unit tests for the methylation-concordance classes.
#
# usage: bash test/unit/concordance/run.sh [scratchDir]
#
# Compiles the main sources and ConcordanceTest.java into <scratchDir> (default: a fresh mktemp
# directory; never target/) and runs every scenario in its own JVM. Exit code 0 = all passed.
#
# Environment:
#   JAVA_HOME     JDK to use (default: Quest openjdk 11)
#   MAIN_CLASSES  directory of already compiled main classes to test instead of compiling the
#                 source tree (used to run the tests against a pre-fix build)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$HERE/../../.." && pwd)
JAVA_HOME=${JAVA_HOME:-/gpfs/software/spack_v20d1/spack/opt/spack/linux-rhel7-x86_64/gcc-10.4.0/openjdk-11.0.17_8-xt4o77czxlhvrycpl4jo3z5swnqihkvz}
SCRATCH=${1:-$(mktemp -d)}
mkdir -p "$SCRATCH/test_classes" "$SCRATCH/work"
DEPS="$R/target/bisulfitehic-0.38-jar-with-dependencies.jar:$R/lib/*:$R/jbwa/jbwa-1.0.0/jbwa.jar"

if [ -z "${MAIN_CLASSES:-}" ]; then
	MAIN_CLASSES="$SCRATCH/main_classes"
	mkdir -p "$MAIN_CLASSES"
	"$JAVA_HOME/bin/javac" -nowarn --release 8 -d "$MAIN_CLASSES" -cp "$DEPS" \
		"$R"/src/main/java/edu/mit/compbio/bisulfitehic/*/*.java 2> "$SCRATCH/javac_main.log" \
		|| { cat "$SCRATCH/javac_main.log"; echo "main sources failed to compile"; exit 2; }
fi

# freshly compiled classes must come before the production jar on the classpath
CP="$MAIN_CLASSES:$DEPS"
"$JAVA_HOME/bin/javac" -nowarn --release 8 -d "$SCRATCH/test_classes" -cp "$CP" "$HERE/ConcordanceTest.java" \
	|| { echo "ConcordanceTest failed to compile"; exit 2; }

status=0
for scenario in helpers random count diffchr short selfpair unpaired global oddbws; do
	"$JAVA_HOME/bin/java" -cp "$SCRATCH/test_classes:$CP" ConcordanceTest "$scenario" "$SCRATCH/work" \
		> "$SCRATCH/$scenario.log" 2>&1
	rc=$?
	grep -E "^(PASS|FAIL|  check failed)" "$SCRATCH/$scenario.log"
	if [ $rc -ne 0 ]; then
		status=1
		grep -q -E "^FAIL" "$SCRATCH/$scenario.log" || { echo "FAIL $scenario (exit $rc, see $SCRATCH/$scenario.log)"; tail -5 "$SCRATCH/$scenario.log"; }
	fi
done
echo "scratch directory: $SCRATCH"
exit $status
