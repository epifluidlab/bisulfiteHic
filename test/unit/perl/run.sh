#!/bin/bash
# Unit test for src/perl/demultiplex_ecker_scWGBS.pl: after removing the 6 bp
# barcode, R1 must keep every remaining base and quality value (an off-by-one used
# to drop the last one), R2 must be untouched, and unmatched barcodes are skipped.
# Usage: bash test/unit/perl/run.sh      Exit status 0 = pass.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
W="$(mktemp -d)"
trap 'rm -rf "$W"' EXIT

printf 'ACGTAC\nTTTTTT\n' > "$W/index.txt"
{
    printf '@r1\nACGTACAAAACCCCGGGGTTTTAC\n+\nJJJJJJABCDEFGHIJKLMNOPQR\n'   # barcode ACGTAC + 18 bases
    printf '@r2\nGGGGGGAAAACCCCGGGGTTTTAC\n+\nJJJJJJABCDEFGHIJKLMNOPQR\n'   # barcode not in the index
} > "$W/in_R1.fq"
{
    printf '@r1\nTTTTGGGGCCCCAAAA\n+\nABCDEFGHIJKLMNOP\n'
    printf '@r2\nTTTTGGGGCCCCAAAA\n+\nABCDEFGHIJKLMNOP\n'
} > "$W/in_R2.fq"

perl "$ROOT/src/perl/demultiplex_ecker_scWGBS.pl" "$W/out" "$W/index.txt" "$W/in_R1.fq" "$W/in_R2.fq" 2> "$W/log" \
    || { echo "FAIL: script exited non-zero"; cat "$W/log"; exit 1; }

fail=0
check() {  # check <description> <expected> <actual>
    if [ "$2" = "$3" ]; then echo "  PASS  $1"; else echo "  FAIL  $1: expected '$2', got '$3'"; fail=1; fi
}
r1="$(zcat "$W/out.ACGTAC.R1.fastq.gz")"
r2="$(zcat "$W/out.ACGTAC.R2.fastq.gz")"
check "R1 sequence keeps all 18 bases after the barcode" "AAAACCCCGGGGTTTTAC" "$(echo "$r1" | sed -n 2p)"
check "R1 quality keeps all 18 values after the barcode" "ABCDEFGHIJKLMNOPQR" "$(echo "$r1" | sed -n 4p)"
check "only the matching read is written to R1" "4" "$(echo "$r1" | wc -l)"
check "R2 sequence is unchanged" "TTTTGGGGCCCCAAAA" "$(echo "$r2" | sed -n 2p)"
check "unused barcode produces an empty file" "0" "$(zcat "$W/out.TTTTTT.R1.fastq.gz" | wc -l)"
if [ $fail -eq 0 ]; then echo "PASS: demultiplex_ecker_scWGBS.pl"; else echo "FAIL: demultiplex_ecker_scWGBS.pl"; fi
exit $fail
