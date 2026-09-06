#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <baseline.bam> <work_dir>" >&2
    exit 1
fi

BASELINE_BAM="$(readlink -f "$1")"
WORK_DIR="$2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPERIMENT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TRUTH_FILE="${EXPERIMENT_DIR}/positive_control_truth.tsv"
REMOVE_SCRIPT="${SCRIPT_DIR}/remove_first_k_mapped.sh"

mkdir -p "$WORK_DIR"
WORK_DIR="$(readlink -f "$WORK_DIR")"

if [[ ! -f "$BASELINE_BAM" ]]; then
    echo "ERROR: baseline BAM not found: $BASELINE_BAM" >&2
    exit 1
fi

if [[ ! -f "${BASELINE_BAM}.bai" ]]; then
    echo "ERROR: matching baseline BAI not found: ${BASELINE_BAM}.bai" >&2
    exit 1
fi

if [[ ! -f "$TRUTH_FILE" ]]; then
    echo "ERROR: truth specification not found: $TRUTH_FILE" >&2
    exit 1
fi

command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools not found" >&2
    exit 1
}

command -v awk >/dev/null 2>&1 || {
    echo "ERROR: awk not found" >&2
    exit 1
}

samtools quickcheck -v "$BASELINE_BAM"

flagstat_mapped() {
    samtools flagstat "$1" \
      | awk '/ mapped \(/ {print $1 + $3; exit}'
}

idxstats_mapped() {
    samtools idxstats "$1" \
      | awk '{mapped += $3} END {print mapped+0}'
}

BASE_NAME="$(basename "$BASELINE_BAM" .bam)"

PC1_BAM="${WORK_DIR}/${BASE_NAME}_PC1_k1000.bam"
PC2_BAM="${WORK_DIR}/${BASE_NAME}_PC2_k10000.bam"

RESULTS="${WORK_DIR}/experiment_2b_results.tsv"

echo "Experiment 2b"
echo "Baseline: $BASELINE_BAM"
echo "Work directory: $WORK_DIR"
echo

echo "===== BASELINE ====="

BASE_F="$(flagstat_mapped "$BASELINE_BAM")"
BASE_I="$(idxstats_mapped "$BASELINE_BAM")"
BASE_DELTA=$((BASE_F - BASE_I))

echo "flagstat mapped = $BASE_F"
echo "idxstats mapped = $BASE_I"
echo "delta = $BASE_DELTA"

echo
echo "===== PC1: k=1000 ====="

"$REMOVE_SCRIPT" "$BASELINE_BAM" 1000 "$PC1_BAM"
samtools quickcheck -v "$PC1_BAM"

PC1_F="$(flagstat_mapped "$PC1_BAM")"

cp "${BASELINE_BAM}.bai" "${PC1_BAM}.bai"

PC1_I="$(idxstats_mapped "$PC1_BAM")"
PC1_DELTA=$((PC1_F - PC1_I))
PC1_EXPECTED=-1000
PC1_RESIDUAL=$((PC1_DELTA - PC1_EXPECTED))

echo "flagstat mapped = $PC1_F"
echo "idxstats mapped = $PC1_I"
echo "delta observed = $PC1_DELTA"
echo "delta expected = $PC1_EXPECTED"
echo "residual = $PC1_RESIDUAL"

echo
echo "===== PC2: k=10000 ====="

"$REMOVE_SCRIPT" "$BASELINE_BAM" 10000 "$PC2_BAM"
samtools quickcheck -v "$PC2_BAM"

PC2_F="$(flagstat_mapped "$PC2_BAM")"

cp "${BASELINE_BAM}.bai" "${PC2_BAM}.bai"

PC2_I="$(idxstats_mapped "$PC2_BAM")"
PC2_DELTA=$((PC2_F - PC2_I))
PC2_EXPECTED=-10000
PC2_RESIDUAL=$((PC2_DELTA - PC2_EXPECTED))

echo "flagstat mapped = $PC2_F"
echo "idxstats mapped = $PC2_I"
echo "delta observed = $PC2_DELTA"
echo "delta expected = $PC2_EXPECTED"
echo "residual = $PC2_RESIDUAL"

cat > "$RESULTS" <<RESULTS_EOF
condition	baseline_mapped_count	removed_mapped_records	flagstat_mapped	idxstats_mapped	delta_expected	delta_observed	residual
baseline	${BASE_F}	0	${BASE_F}	${BASE_I}	0	${BASE_DELTA}	${BASE_DELTA}
PC1_k1000	${BASE_F}	1000	${PC1_F}	${PC1_I}	${PC1_EXPECTED}	${PC1_DELTA}	${PC1_RESIDUAL}
PC2_k10000	${BASE_F}	10000	${PC2_F}	${PC2_I}	${PC2_EXPECTED}	${PC2_DELTA}	${PC2_RESIDUAL}
RESULTS_EOF

echo
echo "===== VALIDATION ====="

EXPECTED_BASE=9061436
EXPECTED_PC1_F=9060436
EXPECTED_PC2_F=9051436

FAILED=0

if [[ "$BASE_F" -ne "$EXPECTED_BASE" || \
      "$BASE_I" -ne "$EXPECTED_BASE" || \
      "$BASE_DELTA" -ne 0 ]]; then
    echo "FAIL: baseline did not reproduce expected concordance." >&2
    FAILED=1
else
    echo "PASS: baseline concordance reproduced."
fi

if [[ "$PC1_F" -ne "$EXPECTED_PC1_F" || \
      "$PC1_I" -ne "$EXPECTED_BASE" || \
      "$PC1_DELTA" -ne -1000 || \
      "$PC1_RESIDUAL" -ne 0 ]]; then
    echo "FAIL: PC1 did not reproduce expected quantitative discrepancy." >&2
    FAILED=1
else
    echo "PASS: PC1 exact discrepancy recovered."
fi

if [[ "$PC2_F" -ne "$EXPECTED_PC2_F" || \
      "$PC2_I" -ne "$EXPECTED_BASE" || \
      "$PC2_DELTA" -ne -10000 || \
      "$PC2_RESIDUAL" -ne 0 ]]; then
    echo "FAIL: PC2 did not reproduce expected quantitative discrepancy." >&2
    FAILED=1
else
    echo "PASS: PC2 exact discrepancy recovered."
fi

if [[ "$FAILED" -ne 0 ]]; then
    echo "Experiment 2b validation FAILED." >&2
    exit 1
fi

echo
echo "Experiment 2b validation PASSED."
echo "Results written to: $RESULTS"
