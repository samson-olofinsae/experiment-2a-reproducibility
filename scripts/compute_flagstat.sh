#!/usr/bin/env bash

set -euo pipefail

# Generate samtools flagstat output for one BAM file used in the
# baseline concordance experiment.
#
# Usage:
#   bash scripts/compute_flagstat.sh input.bam output.txt
#
# Example:
#   bash scripts/compute_flagstat.sh \
#       data/downsampled/SRR10353673_10pct.bam \
#       results/flagstat/SRR10353673_10pct.flagstat.txt

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 input.bam output.txt" >&2
    exit 1
fi

BAM_FILE="$1"
OUTPUT_TXT="$2"

# Confirm that samtools is available.

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools was not found in PATH." >&2
    exit 1
fi

# Confirm that the input BAM exists and is not empty.

if [[ ! -s "${BAM_FILE}" ]]; then
    echo "ERROR: BAM file '${BAM_FILE}' was not found or is empty." >&2
    exit 1
fi

# Ensure that the output directory exists.

mkdir -p "$(dirname "${OUTPUT_TXT}")"

# Generate the raw flagstat report used for mapped-read extraction.

samtools flagstat "${BAM_FILE}" > "${OUTPUT_TXT}"

# Confirm that the output was created and is not empty.

if [[ ! -s "${OUTPUT_TXT}" ]]; then
    echo "ERROR: flagstat output was not created or is empty: ${OUTPUT_TXT}" >&2
    exit 1
fi

echo "Flagstat written to ${OUTPUT_TXT}"