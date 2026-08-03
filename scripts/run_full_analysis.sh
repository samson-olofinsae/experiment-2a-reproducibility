#!/usr/bin/env bash

set -euo pipefail

# Reproduce the complete Experiment 2a analytical workflow.
#
# Starting from the three source BAM files listed in:
#
#   config/experiment_2a_bam_manifest.tsv
#
# this script:
#
#   1. validates the source BAM files using samtools quickcheck;
#   2. verifies or creates BAM indexes;
#   3. independently downsamples each source BAM to 10–90%;
#   4. generates raw samtools flagstat outputs;
#   5. generates raw samtools idxstats outputs;
#   6. creates metrics_table.tsv;
#   7. creates discrepancy_summary.tsv;
#   8. generates manuscript Figures 2–5;
#   9. verifies that all expected outputs were produced.
#
# The script does not download the source BAM files automatically.
# Run scripts/download_bams.sh first, or place the three source BAMs in:
#
#   data/raw_bams/
#
# Usage:
#
#   bash scripts/run_full_analysis.sh
#
# To replace previously generated analytical outputs:
#
#   bash scripts/run_full_analysis.sh --overwrite

MANIFEST="config/experiment_2a_bam_manifest.tsv"

RAW_DIR="data/raw_bams"
DOWNSAMPLED_DIR="data/downsampled"

RESULTS_DIR="results"
FLAGSTAT_DIR="${RESULTS_DIR}/flagstat"
IDXSTATS_DIR="${RESULTS_DIR}/idxstats"
FIGURES_DIR="${RESULTS_DIR}/figures"

METRICS_TABLE="${RESULTS_DIR}/metrics_table.tsv"
SUMMARY_TABLE="${RESULTS_DIR}/discrepancy_summary.tsv"
STRUCTURAL_REPORT="${RESULTS_DIR}/structural_validation.tsv"

VALIDATION_SCRIPT="scripts/validate_pipeline.sh"
DOWNSAMPLE_SCRIPT="scripts/downsample.sh"
FLAGSTAT_SCRIPT="scripts/compute_flagstat.sh"
IDXSTATS_SCRIPT="scripts/compute_idxstats.sh"
TABLE_SCRIPT="scripts/generate_metrics_table.py"
PLOTTING_SCRIPT="scripts/plot_metrics.py"

FRACTIONS=(
    0.1
    0.2
    0.3
    0.4
    0.5
    0.6
    0.7
    0.8
    0.9
)

EXPECTED_SAMPLES=3
EXPECTED_DEPTHS=9
EXPECTED_COMPARISONS=$((EXPECTED_SAMPLES * EXPECTED_DEPTHS))

OVERWRITE=false

# ------------------------------------------------------------------
# Command-line argument handling
# ------------------------------------------------------------------

if [[ $# -gt 1 ]]; then
    echo "Usage: $0 [--overwrite]" >&2
    exit 1
fi

if [[ $# -eq 1 ]]; then
    case "$1" in
        --overwrite)
            OVERWRITE=true
            ;;
        -h|--help)
            echo "Usage: $0 [--overwrite]"
            echo
            echo "Options:"
            echo "  --overwrite   Replace previously generated Experiment 2a outputs."
            echo "  -h, --help    Show this help message."
            exit 0
            ;;
        *)
            echo "ERROR: Unrecognised argument: $1" >&2
            echo "Usage: $0 [--overwrite]" >&2
            exit 1
            ;;
    esac
fi

# ------------------------------------------------------------------
# Confirm execution from the repository root
# ------------------------------------------------------------------

if [[ ! -f "${MANIFEST}" || ! -d "scripts" ]]; then
    echo "ERROR: Required repository files were not found." >&2
    echo "Run this command from the experiment-2a-reproducibility root directory." >&2
    exit 1
fi

# ------------------------------------------------------------------
# Dependency checks
# ------------------------------------------------------------------

for command_name in samtools python3; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "ERROR: Required command was not found in PATH: ${command_name}" >&2
        exit 1
    fi
done

if ! python3 -c "import pandas, matplotlib" >/dev/null 2>&1; then
    echo "ERROR: Required Python libraries pandas and matplotlib are unavailable." >&2
    echo "Create the project environment from environment.yml before continuing." >&2
    exit 1
fi

# ------------------------------------------------------------------
# Required file checks
# ------------------------------------------------------------------

required_files=(
    "${MANIFEST}"
    "${VALIDATION_SCRIPT}"
    "${DOWNSAMPLE_SCRIPT}"
    "${FLAGSTAT_SCRIPT}"
    "${IDXSTATS_SCRIPT}"
    "${TABLE_SCRIPT}"
    "${PLOTTING_SCRIPT}"
)

for required_file in "${required_files[@]}"; do
    if [[ ! -f "${required_file}" ]]; then
        echo "ERROR: Required file not found: ${required_file}" >&2
        exit 1
    fi
done

if [[ ! -d "${RAW_DIR}" ]]; then
    echo "ERROR: Source BAM directory not found: ${RAW_DIR}" >&2
    echo "Run scripts/download_bams.sh or place the Experiment 2a BAM files there." >&2
    exit 1
fi

# ------------------------------------------------------------------
# Read the fixed Experiment 2a sample set
# ------------------------------------------------------------------

mapfile -t SAMPLE_IDS < <(
    tail -n +2 "${MANIFEST}" |
    awk -F $'\t' 'NF > 0 && $1 != "" {print $1}'
)

if [[ ${#SAMPLE_IDS[@]} -ne ${EXPECTED_SAMPLES} ]]; then
    echo "ERROR: Manifest contains ${#SAMPLE_IDS[@]} samples; expected ${EXPECTED_SAMPLES}." >&2
    exit 1
fi

echo "Experiment 2a samples:"

for sample_id in "${SAMPLE_IDS[@]}"; do
    echo "  - ${sample_id}"

    source_bam="${RAW_DIR}/${sample_id}.bam"

    if [[ ! -s "${source_bam}" ]]; then
        echo "ERROR: Required source BAM is missing or empty: ${source_bam}" >&2
        exit 1
    fi
done

# ------------------------------------------------------------------
# Protect against stale analytical outputs
# ------------------------------------------------------------------

generated_paths=(
    "${DOWNSAMPLED_DIR}"
    "${FLAGSTAT_DIR}"
    "${IDXSTATS_DIR}"
    "${FIGURES_DIR}"
    "${METRICS_TABLE}"
    "${SUMMARY_TABLE}"
)

existing_outputs=false

for generated_path in "${generated_paths[@]}"; do
    if [[ -e "${generated_path}" ]]; then
        existing_outputs=true
        break
    fi
done

if [[ "${existing_outputs}" == true && "${OVERWRITE}" == false ]]; then
    echo "ERROR: Previously generated analytical outputs already exist." >&2
    echo "To replace them, run:" >&2
    echo >&2
    echo "  bash scripts/run_full_analysis.sh --overwrite" >&2
    exit 1
fi

if [[ "${OVERWRITE}" == true ]]; then
    echo "Removing previously generated analytical outputs..."

    rm -rf \
        "${DOWNSAMPLED_DIR}" \
        "${FLAGSTAT_DIR}" \
        "${IDXSTATS_DIR}" \
        "${FIGURES_DIR}"

    rm -f \
        "${METRICS_TABLE}" \
        "${SUMMARY_TABLE}"
fi

mkdir -p \
    "${DOWNSAMPLED_DIR}" \
    "${FLAGSTAT_DIR}" \
    "${IDXSTATS_DIR}" \
    "${FIGURES_DIR}"

# ------------------------------------------------------------------
# Step 1: structural validation of source BAM files
# ------------------------------------------------------------------

echo
echo "============================================================"
echo "Step 1 of 5: Validating source BAM files"
echo "============================================================"

bash "${VALIDATION_SCRIPT}"

if [[ ! -s "${STRUCTURAL_REPORT}" ]]; then
    echo "ERROR: Structural-validation report was not produced." >&2
    exit 1
fi

validation_pass_count="$(
    awk -F $'\t' '
        NR > 1 && $5 == "PASS" {
            count++
        }
        END {
            print count + 0
        }
    ' "${STRUCTURAL_REPORT}"
)"

if [[ "${validation_pass_count}" -ne "${EXPECTED_SAMPLES}" ]]; then
    echo "ERROR: Only ${validation_pass_count} source BAM files passed validation;" >&2
    echo "expected ${EXPECTED_SAMPLES}." >&2
    exit 1
fi

# ------------------------------------------------------------------
# Steps 2 and 3: downsampling and raw metric extraction
# ------------------------------------------------------------------

echo
echo "============================================================"
echo "Steps 2–3 of 5: Downsampling and metric extraction"
echo "============================================================"

comparison_count=0

for sample_id in "${SAMPLE_IDS[@]}"; do
    source_bam="${RAW_DIR}/${sample_id}.bam"

    echo
    echo "Processing source BAM: ${sample_id}"

    for fraction in "${FRACTIONS[@]}"; do
        depth_percent="$(
            python3 - "${fraction}" <<'PY'
import sys

fraction = float(sys.argv[1])
print(int(round(fraction * 100)))
PY
        )"

        output_stem="${sample_id}_${depth_percent}pct"

        downsampled_bam="${DOWNSAMPLED_DIR}/${output_stem}.bam"
        flagstat_output="${FLAGSTAT_DIR}/${output_stem}.flagstat.txt"
        idxstats_output="${IDXSTATS_DIR}/${output_stem}.idxstats.txt"

        echo
        echo "  Depth: ${depth_percent}%"

        bash \
            "${DOWNSAMPLE_SCRIPT}" \
            "${source_bam}" \
            "${fraction}" \
            "${downsampled_bam}"

        bash \
            "${FLAGSTAT_SCRIPT}" \
            "${downsampled_bam}" \
            "${flagstat_output}"

        bash \
            "${IDXSTATS_SCRIPT}" \
            "${downsampled_bam}" \
            "${idxstats_output}"

        comparison_count=$((comparison_count + 1))
    done
done

if [[ "${comparison_count}" -ne "${EXPECTED_COMPARISONS}" ]]; then
    echo "ERROR: Generated ${comparison_count} comparisons;" >&2
    echo "expected ${EXPECTED_COMPARISONS}." >&2
    exit 1
fi

downsampled_count="$(
    find "${DOWNSAMPLED_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.bam' |
    wc -l
)"

flagstat_count="$(
    find "${FLAGSTAT_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.flagstat.txt' |
    wc -l
)"

idxstats_count="$(
    find "${IDXSTATS_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.idxstats.txt' |
    wc -l
)"

if [[ "${downsampled_count}" -ne "${EXPECTED_COMPARISONS}" ]]; then
    echo "ERROR: Found ${downsampled_count} downsampled BAM files;" >&2
    echo "expected ${EXPECTED_COMPARISONS}." >&2
    exit 1
fi

if [[ "${flagstat_count}" -ne "${EXPECTED_COMPARISONS}" ]]; then
    echo "ERROR: Found ${flagstat_count} flagstat outputs;" >&2
    echo "expected ${EXPECTED_COMPARISONS}." >&2
    exit 1
fi

if [[ "${idxstats_count}" -ne "${EXPECTED_COMPARISONS}" ]]; then
    echo "ERROR: Found ${idxstats_count} idxstats outputs;" >&2
    echo "expected ${EXPECTED_COMPARISONS}." >&2
    exit 1
fi

# ------------------------------------------------------------------
# Step 4: table generation
# ------------------------------------------------------------------

echo
echo "============================================================"
echo "Step 4 of 5: Generating benchmarking tables"
echo "============================================================"

python3 "${TABLE_SCRIPT}"

if [[ ! -s "${METRICS_TABLE}" ]]; then
    echo "ERROR: Metrics table was not generated: ${METRICS_TABLE}" >&2
    exit 1
fi

if [[ ! -s "${SUMMARY_TABLE}" ]]; then
    echo "ERROR: Summary table was not generated: ${SUMMARY_TABLE}" >&2
    exit 1
fi

metrics_row_count="$(
    awk 'END {print NR - 1}' "${METRICS_TABLE}"
)"

summary_row_count="$(
    awk 'END {print NR - 1}' "${SUMMARY_TABLE}"
)"

if [[ "${metrics_row_count}" -ne "${EXPECTED_COMPARISONS}" ]]; then
    echo "ERROR: Metrics table contains ${metrics_row_count} data rows;" >&2
    echo "expected ${EXPECTED_COMPARISONS}." >&2
    exit 1
fi

if [[ "${summary_row_count}" -ne "${EXPECTED_DEPTHS}" ]]; then
    echo "ERROR: Summary table contains ${summary_row_count} data rows;" >&2
    echo "expected ${EXPECTED_DEPTHS}." >&2
    exit 1
fi

# ------------------------------------------------------------------
# Step 5: figure generation
# ------------------------------------------------------------------

echo
echo "============================================================"
echo "Step 5 of 5: Generating publication figures"
echo "============================================================"

python3 "${PLOTTING_SCRIPT}"

expected_figure_outputs=(
    "${FIGURES_DIR}/figure2_mean_mapped_reads_by_depth.png"
    "${FIGURES_DIR}/figure2_mean_mapped_reads_by_depth.pdf"
    "${FIGURES_DIR}/figure3_flagstat_vs_idxstats_agreement.png"
    "${FIGURES_DIR}/figure3_flagstat_vs_idxstats_agreement.pdf"
    "${FIGURES_DIR}/figure4_mapped_difference_pct_by_depth.png"
    "${FIGURES_DIR}/figure4_mapped_difference_pct_by_depth.pdf"
    "${FIGURES_DIR}/figure5_sample_mapped_reads_by_depth.png"
    "${FIGURES_DIR}/figure5_sample_mapped_reads_by_depth.pdf"
    "${FIGURES_DIR}/figure_summary_stats.txt"
)

for output_file in "${expected_figure_outputs[@]}"; do
    if [[ ! -s "${output_file}" ]]; then
        echo "ERROR: Expected figure output is missing or empty: ${output_file}" >&2
        exit 1
    fi
done

# ------------------------------------------------------------------
# Completion summary
# ------------------------------------------------------------------

echo
echo "============================================================"
echo "Experiment 2a reproduction completed successfully"
echo "============================================================"
echo
echo "Source BAM files validated: ${validation_pass_count}"
echo "Downsampled BAM files generated: ${downsampled_count}"
echo "flagstat outputs generated: ${flagstat_count}"
echo "idxstats outputs generated: ${idxstats_count}"
echo "Benchmarking comparisons: ${metrics_row_count}"
echo "Depth-level summary rows: ${summary_row_count}"
echo
echo "Primary outputs:"
echo "  ${METRICS_TABLE}"
echo "  ${SUMMARY_TABLE}"
echo "  ${FIGURES_DIR}/"
echo
echo "Canonical outputs for comparison are available in:"
echo "  expected_results/"