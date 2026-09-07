#!/usr/bin/env bash

set -euo pipefail

# Structural validation of the source BAM files used in the
# baseline concordance experiment.
#
# This script:
#   1. Reads the fixed source-BAM manifest.
#   2. Confirms that each expected BAM file is present.
#   3. Runs samtools quickcheck to assess BAM structure.
#   4. Confirms that a compatible BAM index is present or creates one.
#   5. Moves structurally invalid BAM files and associated indexes
#      into the quarantine directory.
#   6. Writes tab-separated structural-validation reports.
#
# Important:
#   This script verifies BAM structure and index availability.
#   It does not determine whether an existing BAM index represents
#   the current contents of the BAM. Stale-index behaviour is evaluated
#   separately in the prespecified positive-control experiment.
#
# Expected inputs:
#   config/experiment_2a_bam_manifest.tsv
#   data/raw_bams/<sample_id>.bam
#
# Outputs:
#   results/structural_validation.tsv
#   results/structural_validation_failures.tsv
#   data/quarantine/

MANIFEST="config/experiment_2a_bam_manifest.tsv"

RAW_DIR="data/raw_bams"
QUARANTINE_DIR="data/quarantine"
RESULTS_DIR="results"

VALIDATION_REPORT="${RESULTS_DIR}/structural_validation.tsv"
FAILURE_REPORT="${RESULTS_DIR}/structural_validation_failures.tsv"

# ------------------------------------------------------------------
# Dependency and input checks
# ------------------------------------------------------------------

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools was not found in PATH." >&2
    exit 1
fi

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}" >&2
    exit 1
fi

if [[ ! -d "${RAW_DIR}" ]]; then
    echo "ERROR: Raw BAM directory not found: ${RAW_DIR}" >&2
    echo "Run scripts/download_bams.sh first or place the validated source BAMs there." >&2
    exit 1
fi

mkdir -p "${QUARANTINE_DIR}" "${RESULTS_DIR}"

# ------------------------------------------------------------------
# Initialise reports
# ------------------------------------------------------------------

printf \
    "sample_id\tbam_file\tquickcheck_status\tindex_status\toverall_status\tnotes\n" \
    > "${VALIDATION_REPORT}"

printf \
    "sample_id\tbam_file\tfailure_type\taction_taken\tnotes\n" \
    > "${FAILURE_REPORT}"

echo "Reading manifest: ${MANIFEST}"
echo "Validating source BAM files in: ${RAW_DIR}"

validated_count=0
failed_count=0

# ------------------------------------------------------------------
# Validate each source BAM
# ------------------------------------------------------------------

while IFS=$'\t' read -r \
    sample_id \
    study_accession \
    cancer_type \
    assay_type \
    source \
    file_url \
    expected_md5 \
    notes
do
    [[ "${sample_id}" == "sample_id" ]] && continue
    [[ -z "${sample_id:-}" ]] && continue

    bam_file="${RAW_DIR}/${sample_id}.bam"

    # samtools commonly recognises either:
    #   sample.bam.bai
    # or:
    #   sample.bai
    primary_index="${bam_file}.bai"
    alternate_index="${bam_file%.bam}.bai"

    echo "----------------------------------------"
    echo "Sample: ${sample_id}"
    echo "BAM: ${bam_file}"

    # --------------------------------------------------------------
    # Confirm that the BAM exists and is not empty
    # --------------------------------------------------------------

    if [[ ! -s "${bam_file}" ]]; then
        echo "FAIL: BAM file is missing or empty."

        printf \
            "%s\t%s\tNOT_RUN\tNOT_RUN\tFAIL\tBAM file missing or empty\n" \
            "${sample_id}" \
            "${bam_file}" \
            >> "${VALIDATION_REPORT}"

        printf \
            "%s\t%s\tMISSING_OR_EMPTY\tNo file moved\tExpected BAM was absent or empty\n" \
            "${sample_id}" \
            "${bam_file}" \
            >> "${FAILURE_REPORT}"

        failed_count=$((failed_count + 1))
        continue
    fi

    # --------------------------------------------------------------
    # Structural validation with samtools quickcheck
    # --------------------------------------------------------------

    quickcheck_output=""

    if quickcheck_output="$(samtools quickcheck -v "${bam_file}" 2>&1)"; then
        quickcheck_status="PASS"
        echo "samtools quickcheck: PASS"
    else
        quickcheck_status="FAIL"

        if [[ -z "${quickcheck_output}" ]]; then
            quickcheck_output="samtools quickcheck returned a non-zero exit status"
        fi

        echo "samtools quickcheck: FAIL"
        echo "${quickcheck_output}" >&2

        timestamp="$(date +%Y%m%d_%H%M%S)"
        quarantine_bam="${QUARANTINE_DIR}/${sample_id}_quickcheck_failed_${timestamp}.bam"

        mv "${bam_file}" "${quarantine_bam}"

        moved_indexes="none"

        if [[ -e "${primary_index}" ]]; then
            quarantine_index="${quarantine_bam}.bai"
            mv "${primary_index}" "${quarantine_index}"
            moved_indexes="${quarantine_index}"
        fi

        if [[ "${alternate_index}" != "${primary_index}" && -e "${alternate_index}" ]]; then
            alternate_quarantine_index="${QUARANTINE_DIR}/${sample_id}_quickcheck_failed_${timestamp}.bai"
            mv "${alternate_index}" "${alternate_quarantine_index}"

            if [[ "${moved_indexes}" == "none" ]]; then
                moved_indexes="${alternate_quarantine_index}"
            else
                moved_indexes="${moved_indexes};${alternate_quarantine_index}"
            fi
        fi

        printf \
            "%s\t%s\tFAIL\tNOT_RUN\tFAIL\tStructural validation failed; BAM moved to %s\n" \
            "${sample_id}" \
            "${bam_file}" \
            "${quarantine_bam}" \
            >> "${VALIDATION_REPORT}"

        printf \
            "%s\t%s\tQUICKCHECK_FAILED\tMoved BAM to %s\t%s; indexes moved: %s\n" \
            "${sample_id}" \
            "${bam_file}" \
            "${quarantine_bam}" \
            "${quickcheck_output//$'\t'/ }" \
            "${moved_indexes}" \
            >> "${FAILURE_REPORT}"

        failed_count=$((failed_count + 1))
        continue
    fi

    # --------------------------------------------------------------
    # BAM index verification or creation
    # --------------------------------------------------------------

    if [[ -s "${primary_index}" ]]; then
        index_file="${primary_index}"
        index_status="PRESENT"
        echo "BAM index: PRESENT (${index_file})"

    elif [[ -s "${alternate_index}" ]]; then
        index_file="${alternate_index}"
        index_status="PRESENT"
        echo "BAM index: PRESENT (${index_file})"

    else
        echo "BAM index: NOT FOUND"
        echo "Creating BAM index with samtools index..."

        if samtools index "${bam_file}"; then
            if [[ -s "${primary_index}" ]]; then
                index_file="${primary_index}"
                index_status="CREATED"
            elif [[ -s "${alternate_index}" ]]; then
                index_file="${alternate_index}"
                index_status="CREATED"
            else
                index_file="NA"
                index_status="FAIL"
            fi
        else
            index_file="NA"
            index_status="FAIL"
        fi

        if [[ "${index_status}" == "CREATED" ]]; then
            echo "BAM index: CREATED (${index_file})"
        else
            echo "BAM index: FAIL" >&2
        fi
    fi

    # --------------------------------------------------------------
    # Record final outcome
    # --------------------------------------------------------------

    if [[ "${quickcheck_status}" == "PASS" && \
          ( "${index_status}" == "PRESENT" || "${index_status}" == "CREATED" ) ]]; then

        overall_status="PASS"
        validation_notes="BAM passed samtools quickcheck; index ${index_status,,}: ${index_file}"

        printf \
            "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${sample_id}" \
            "${bam_file}" \
            "${quickcheck_status}" \
            "${index_status}" \
            "${overall_status}" \
            "${validation_notes}" \
            >> "${VALIDATION_REPORT}"

        validated_count=$((validated_count + 1))
    else
        overall_status="FAIL"
        validation_notes="BAM passed quickcheck but index verification or creation failed"

        printf \
            "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${sample_id}" \
            "${bam_file}" \
            "${quickcheck_status}" \
            "${index_status}" \
            "${overall_status}" \
            "${validation_notes}" \
            >> "${VALIDATION_REPORT}"

        printf \
            "%s\t%s\tINDEX_FAILED\tBAM retained; index could not be verified or created\t%s\n" \
            "${sample_id}" \
            "${bam_file}" \
            "${validation_notes}" \
            >> "${FAILURE_REPORT}"

        failed_count=$((failed_count + 1))
    fi

done < "${MANIFEST}"

# ------------------------------------------------------------------
# Final summary and exit status
# ------------------------------------------------------------------

echo "----------------------------------------"
echo "Structural validation complete."
echo "Source BAM files passed: ${validated_count}"
echo "Source BAM files failed: ${failed_count}"
echo "Validation report: ${VALIDATION_REPORT}"
echo "Failure report: ${FAILURE_REPORT}"
echo "Quarantine directory: ${QUARANTINE_DIR}"

if [[ "${failed_count}" -gt 0 ]]; then
    echo "ERROR: One or more baseline source BAM files failed structural validation." >&2
    exit 1
fi

echo "All baseline source BAM files passed structural validation."