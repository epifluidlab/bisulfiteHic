#!/bin/bash
# End-to-end regression tests for bisulfitehic (Bhmem aligner + helper scripts).
#
# Builds a tiny synthetic bisulfite reference, simulates directional / PBAT Hi-C
# read pairs, runs Bhmem from the built jar and checks the BAM. Takes ~1 minute.
#
# Usage:   bash test/regression/run_tests.sh
# Env overrides:
#   JAR            jar to test (default: target/bisulfitehic-0.38-jar-with-dependencies.jar)
#   WORK           scratch directory (default: test/regression/work, git-ignored)
#   JAVA_MODULE    module providing java if java is not on PATH (default: java/jdk-17.0.2+8)
#   PYSAM_PYTHON   python with pysam for the mh_hic_filter test
# Requires: java, samtools, python3, and the bundled bwa + libbwajni.so that
# install.sh builds under jbwa/jbwa-1.0.0/.
# Exit status: 0 if every test passes, 1 otherwise.

set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
JAR="${JAR:-$ROOT/target/bisulfitehic-0.38-jar-with-dependencies.jar}"
WORK="${WORK:-$HERE/work}"
JAVA_MODULE="${JAVA_MODULE:-java/jdk-17.0.2+8}"
JBWA="$ROOT/jbwa/jbwa-1.0.0"
BWA="$(ls "$JBWA"/bwa-*/bwa 2>/dev/null | head -n1)"
PYSAM_PYTHON="${PYSAM_PYTHON:-/projects/b1198/epifluidlab/yoshii/software/conda/envs/scnomehic/bin/python}"

if ! command -v java >/dev/null 2>&1; then
    [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
    module load "$JAVA_MODULE" || { echo "cannot find java (tried module $JAVA_MODULE)"; exit 1; }
fi
for tool in java samtools python3; do
    command -v "$tool" >/dev/null 2>&1 || { echo "missing required tool: $tool"; exit 1; }
done
[ -x "$BWA" ] || { echo "bundled bwa not found under $JBWA (run install.sh first)"; exit 1; }
[ -f "$JBWA/src/main/native/libbwajni.so" ] || { echo "libbwajni.so not built (run install.sh first)"; exit 1; }
[ -f "$JAR" ] || { echo "jar not found: $JAR"; exit 1; }

echo "jar : $JAR"
echo "work: $WORK"
rm -rf "$WORK"
mkdir -p "$WORK/out"

# ── test data ────────────────────────────────────────────────────────────────
python3 "$HERE/make_test_data.py" "$WORK" || exit 1
REF="$WORK/ref/genome.fa"
samtools faidx "$REF" || exit 1
samtools dict "$REF" -o "$WORK/ref/genome.dict" 2>/dev/null || exit 1
for conv in CT GA; do
    "$BWA" index "$WORK/ref/Bisulfite_Genome/${conv}_conversion/genome_mfa.${conv}_conversion.fa" \
        > "$WORK/out/bwa_index_${conv}.log" 2>&1 || { echo "bwa index failed for $conv"; exit 1; }
done

bhmem() {  # bhmem <ref.fa> <out.bam> <r1> <r2> <log> [extra Bhmem options...]
    local ref="$1" out="$2" r1="$3" r2="$4" log="$5"
    shift 5
    java -Xmx2G -Djava.library.path="$JBWA/src/main/native" \
        -cp "$JAR:$JBWA/jbwa.jar" \
        main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem \
        "$ref" "$out" "$r1" "$r2" -t 2 -rgId test -rgSm test -buffer 1000 "$@" > "$log" 2>&1
}

FAILED=0
report() {  # report <status> <test name>
    if [ "$1" -eq 0 ]; then echo "PASS: $2"; else echo "FAIL: $2"; FAILED=1; fi
    echo
}

READS="$WORK/reads"
OUT="$WORK/out"

# ── 1. directional library ───────────────────────────────────────────────────
echo "== test 1: directional library, default mode"
bhmem "$REF" "$OUT/dir.bam" "$READS/dir_R1.fq" "$READS/dir_R2.fq" "$OUT/dir.log" -outputMateDiffChr
samtools view "$OUT/dir.bam" 2>/dev/null \
    | python3 "$HERE/check_bam.py" --ref "$REF" --r1 "$READS/dir_R1.fq" --r2 "$READS/dir_R2.fq" --min-pairs 400
report $? "directional BAM has consistent SEQ, QUAL and strand"

# ── 2. PBAT non-directional library (production scNOMe-HiC flags) ────────────
echo "== test 2: PBAT library, -pbat -nonDirectional"
bhmem "$REF" "$OUT/pbat.bam" "$READS/pbat_R1.fq" "$READS/pbat_R2.fq" "$OUT/pbat.log" \
    -pbat -nonDirectional -outputMateDiffChr
samtools view "$OUT/pbat.bam" 2>/dev/null \
    | python3 "$HERE/check_bam.py" --ref "$REF" --r1 "$READS/pbat_R1.fq" --r2 "$READS/pbat_R2.fq" --pbat --min-pairs 400
report $? "PBAT BAM has consistent SEQ, QUAL and strand"

# ── 3. -enzymeList decides between MAPQ-tied candidate pairs ─────────────────
echo "== test 3: -enzymeList"
for label in A B; do
    bhmem "$REF" "$OUT/enz_$label.bam" "$READS/enz_R1.fq" "$READS/enz_R2.fq" "$OUT/enz_$label.log" \
        -pbat -nonDirectional -outputMateDiffChr -enzymeList "$READS/enz_$label.bed"
    samtools view "$OUT/enz_$label.bam" > "$OUT/enz_$label.sam" 2>/dev/null
done
python3 "$HERE/check_enzyme.py" "$READS/enz_truth.tsv" "$OUT/enz_A.sam" "$OUT/enz_B.sam"
report $? "-enzymeList selects the pair supported by restriction sites"

# ── 4. sequence dictionary lookup ────────────────────────────────────────────
echo "== test 4: sequence dictionary lookup"
# 4a: a directory name matching the regex '.fa' ("alfa") must not corrupt the .dict path.
mkdir -p "$WORK/alfa"
cp -r "$WORK/ref/." "$WORK/alfa/"
bhmem "$WORK/alfa/genome.fa" "$OUT/dict_a.bam" "$READS/enz_R1.fq" "$READS/enz_R2.fq" "$OUT/dict_a.log" -pbat -outputMateDiffChr
n=$(samtools view -c "$OUT/dict_a.bam" 2>/dev/null || echo 0)
[ "${n:-0}" -gt 0 ]; s=$?
[ $s -eq 0 ] || { echo "  no alignments written; see $OUT/dict_a.log"; tail -n 3 "$OUT/dict_a.log" | sed 's/^/    /'; }
report $s "reference under a directory called 'alfa' is handled"
# 4b: without a .dict the dictionary must be built from the .fai.
mkdir -p "$WORK/nodict"
cp -r "$WORK/ref/." "$WORK/nodict/"
rm -f "$WORK/nodict/genome.dict"
bhmem "$WORK/nodict/genome.fa" "$OUT/dict_b.bam" "$READS/enz_R1.fq" "$READS/enz_R2.fq" "$OUT/dict_b.log" -pbat -outputMateDiffChr
n=$(samtools view -c "$OUT/dict_b.bam" 2>/dev/null || echo 0)
[ "${n:-0}" -gt 0 ]; s=$?
[ $s -eq 0 ] || { echo "  no alignments written; see $OUT/dict_b.log"; tail -n 3 "$OUT/dict_b.log" | sed 's/^/    /'; }
report $s "missing .dict falls back to the .fai"
# 4c: reference given as a bare file name (no directory component).
( cd "$WORK/ref" && bhmem genome.fa "$OUT/dict_c.bam" "$READS/enz_R1.fq" "$READS/enz_R2.fq" "$OUT/dict_c.log" -pbat -outputMateDiffChr )
n=$(samtools view -c "$OUT/dict_c.bam" 2>/dev/null || echo 0)
[ "${n:-0}" -gt 0 ]; s=$?
[ $s -eq 0 ] || { echo "  no alignments written; see $OUT/dict_c.log"; tail -n 3 "$OUT/dict_c.log" | sed 's/^/    /'; }
report $s "reference given as a bare file name is handled"

# ── 5. mh_hic_filter.py counts long-range cis pairs whose TLEN is 0 ──────────
echo "== test 5: mh_hic_filter.py"
if "$PYSAM_PYTHON" -c "import pysam" >/dev/null 2>&1; then
    {
        printf '@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:2000000\n'
        # three cis pairs 50 kb apart; Bhmem writes TLEN=0 for such chimeric (RF) pairs
        for i in 1 2 3; do
            p1=$((100000 * i)); p2=$((p1 - 50000))
            printf 'p%d\t67\tchr1\t%d\t60\t10M\t=\t%d\t0\tACGTACGTAC\tJJJJJJJJJJ\n' "$i" "$p1" "$p2"
            printf 'p%d\t147\tchr1\t%d\t60\t10M\t=\t%d\t0\tACGTACGTAC\tJJJJJJJJJJ\n' "$i" "$p2" "$p1"
        done
    } > "$OUT/filter.sam"
    samtools view -b "$OUT/filter.sam" -o "$OUT/filter.bam" 2>/dev/null
    "$PYSAM_PYTHON" "$ROOT/src/python/mh_hic_filter.py" --in_cram "$OUT/filter.bam" --out_summary "$OUT/filter.txt" 2> "$OUT/filter.log"
    sed 's/^/  /' "$OUT/filter.txt" 2>/dev/null
    grep -q "^UniqMappedMapQ30NoPcrCisMore10kb:[[:space:]]3$" "$OUT/filter.txt" 2>/dev/null \
        && grep -q "^PerChromCisLongFilter:[[:space:]]0$" "$OUT/filter.txt"
    report $? "mh_hic_filter.py measures cis distance from positions, not TLEN"
else
    echo "SKIPPED: no python with pysam at $PYSAM_PYTHON (set PYSAM_PYTHON)"
    FAILED=1
fi

if [ "$FAILED" -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "SOME TESTS FAILED"; fi
exit "$FAILED"
