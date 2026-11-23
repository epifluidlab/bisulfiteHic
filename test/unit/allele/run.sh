#!/bin/bash
# Compile and run the allele-tool unit tests in this folder against the CURRENT sources
# under src/main/java (not against the classes inside the production jar).
#
#   bash test/unit/allele/run.sh [scratch_dir]
#
# Nothing is written into the repository: classes and synthetic data go to scratch_dir
# (default: a new mktemp directory). Exit status 0 = all tests passed, 1 = a test failed,
# 2 = compilation failed.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R="$(cd "$HERE/../../.." && pwd)"
: "${JAVA_HOME:=/gpfs/software/spack_v20d1/spack/opt/spack/linux-rhel7-x86_64/gcc-10.4.0/openjdk-11.0.17_8-xt4o77czxlhvrycpl4jo3z5swnqihkvz}"
OUT="${1:-$(mktemp -d)}"
mkdir -p "$OUT/classes" "$OUT/data"

JAR="$(ls "$R"/target/bisulfitehic-*-jar-with-dependencies.jar | head -n 1)"
CP="$JAR:$R/lib/*:$R/jbwa/jbwa-1.0.0/jbwa.jar"

# -sourcepath + -Xprefer:source: the tool classes are compiled from src/, never taken from the jar
"$JAVA_HOME/bin/javac" -nowarn -proc:none --release 8 -Xprefer:source -sourcepath "$R/src" \
	-d "$OUT/classes" -cp "$CP" "$HERE"/*.java || { echo "COMPILE FAILED"; exit 2; }

rc=0
for t in LongRangeAsmTest AlleImbalanceExprTest SepAlleleMethySCTest SepAlleleMethySC2Test EnhanceHicMatrixTest; do
	echo "=== $t ==="
	"$JAVA_HOME/bin/java" -cp "$OUT/classes:$CP" "$t" "$OUT/data/$t" 2> "$OUT/data/$t.stderr.log" || rc=1
done
echo "scratch directory: $OUT (tool stderr in data/*.stderr.log)"
exit $rc
