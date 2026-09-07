#!/usr/bin/env bash

set -euo pipefail

# Generate samtools idxstats output for one BAM file used in the
# baseline concordance experiment.
#
# Usage:
#   bash scripts/compute_idxstats.sh input.bam output.txt
#
# Example:
#   bash scripts/compute_idxstats.sh \
#       data/downsampled/SRR10353673_10pct.bam \
#       results/idxstats/SRR10353673_10pct.idxstats.txt

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

PRIMARY_INDEX="${BAM_FILE}.bai"
ALTERNATE_INDEX="${BAM_FILE%.bam}.bai"

# Ensure that a BAM index is available.

if [[ -s "${PRIMARY_INDEX}" ]]; then
    INDEX_FILE="${PRIMARY_INDEX}"
elif [[ -s "${ALTERNATE_INDEX}" ]]; then
    INDEX_FILE="${ALTERNATE_INDEX}"
else
    echo "BAM index not found for ${BAM_FILE}. Creating index..."

    samtools index "${BAM_FILE}"

    if [[ -s "${PRIMARY_INDEX}" ]]; then
        INDEX_FILE="${PRIMARY_INDEX}"
    elif [[ -s "${ALTERNATE_INDEX}" ]]; then
        INDEX_FILE="${ALTERNATE_INDEX}"
    else
        echo "ERROR: BAM index could not be created for ${BAM_FILE}." >&2
        exit 1
    fi
fi

# Generate the raw idxstats report used for mapped-read extraction.

samtools idxstats "${BAM_FILE}" > "${OUTPUT_TXT}"

# Confirm that the output was created and is not empty.

if [[ ! -s "${OUTPUT_TXT}" ]]; then
    echo "ERROR: idxstats output was not created or is empty: ${OUTPUT_TXT}" >&2
    exit 1
fi

echo "Idxstats written to ${OUTPUT_TXT}"
echo "BAM index used: ${INDEX_FILE}"