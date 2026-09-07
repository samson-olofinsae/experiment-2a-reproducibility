#!/usr/bin/env bash

set -euo pipefail

# Generate one deterministic downsampled BAM for the baseline
# concordance experiment.
#
# Usage:
#   bash scripts/downsample.sh input.bam fraction output.bam
#
# Example:
#   bash scripts/downsample.sh sample.bam 0.1 sample_10pct.bam
#
# Reproducibility note:
#
# samtools interprets the -s argument using the form INT.FRAC,
# where INT specifies the random seed and FRAC specifies the
# subsampling fraction.
#
# The baseline concordance experiment uses fractions 0.1 through
# 0.9. The integer component is therefore 0 for every downsampling
# operation, giving a fixed seed of 0 and preserving reproducibility
# of the canonical baseline outputs.
#
# Each fraction is generated directly from the original source BAM,
# rather than sequentially from another downsampled BAM. With the
# fixed seed, the resulting fractions form reproducible nested
# alignment sets across increasing depth conditions.

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 input.bam fraction output.bam" >&2
    exit 1
fi

INPUT_BAM="$1"
FRACTION="$2"
OUTPUT_BAM="$3"

if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "ERROR: Input BAM '${INPUT_BAM}' not found." >&2
    exit 1
fi

# Confirm that the requested fraction is a numeric value greater
# than 0 and less than 1.

python3 - "${FRACTION}" <<'PY'
import sys

try:
    fraction = float(sys.argv[1])
except ValueError:
    print(
        f"ERROR: Downsampling fraction must be numeric: {sys.argv[1]}",
        file=sys.stderr,
    )
    raise SystemExit(1)

if not 0.0 < fraction < 1.0:
    print(
        f"ERROR: Downsampling fraction must be greater than 0 and less than 1: "
        f"{fraction}",
        file=sys.stderr,
    )
    raise SystemExit(1)
PY

# Confirm that samtools is available.

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools was not found in PATH." >&2
    exit 1
fi

# Ensure that the output directory exists.

mkdir -p "$(dirname "${OUTPUT_BAM}")"

# Downsample using the fixed-seed configuration used to generate
# the canonical baseline concordance outputs.

samtools view \
    -@ 2 \
    -s "${FRACTION}" \
    -b "${INPUT_BAM}" \
    -o "${OUTPUT_BAM}"

# Confirm that the output BAM was created and is not empty.

if [[ ! -s "${OUTPUT_BAM}" ]]; then
    echo "ERROR: Downsampled BAM was not created or is empty: ${OUTPUT_BAM}" >&2
    exit 1
fi

# Index the downsampled BAM.

samtools index "${OUTPUT_BAM}"

# Confirm that the BAM index was created.

if [[ ! -s "${OUTPUT_BAM}.bai" ]]; then
    echo "ERROR: BAM index was not created: ${OUTPUT_BAM}.bai" >&2
    exit 1
fi

echo \
    "Downsampled ${INPUT_BAM} -> ${OUTPUT_BAM} at fraction ${FRACTION} " \
    "(samtools seed 0)."