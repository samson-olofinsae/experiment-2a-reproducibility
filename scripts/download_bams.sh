#!/usr/bin/env bash

set -euo pipefail

# Download and MD5-validate the three source BAM files used in the
# baseline concordance experiment.
#
# Input manifest:
#   config/experiment_2a_bam_manifest.tsv
#
# Outputs:
#   data/raw_bams/
#   data/quarantine/
#   results/download_logs/
#   results/download_failures.tsv
#   results/download_validation.tsv

MANIFEST="config/experiment_2a_bam_manifest.tsv"

OUTDIR="data/raw_bams"
QUARANTINE_DIR="data/quarantine"
LOGDIR="results/download_logs"

FAIL_LOG="results/download_failures.tsv"
VALIDATION_LOG="results/download_validation.tsv"

# ------------------------------------------------------------------
# Dependency and input checks
# ------------------------------------------------------------------

if ! command -v aria2c >/dev/null 2>&1; then
    echo "ERROR: aria2c was not found in PATH." >&2
    exit 1
fi

if ! command -v md5sum >/dev/null 2>&1; then
    echo "ERROR: md5sum was not found in PATH." >&2
    exit 1
fi

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}" >&2
    exit 1
fi

mkdir -p \
    "${OUTDIR}" \
    "${QUARANTINE_DIR}" \
    "${LOGDIR}" \
    results

# ------------------------------------------------------------------
# Initialise output logs
# ------------------------------------------------------------------

printf \
    "sample_id\tstudy_accession\tfile_url\texpected_md5\tobserved_md5\tstatus\tnotes\n" \
    > "${FAIL_LOG}"

printf \
    "sample_id\tstudy_accession\tfile_url\texpected_md5\tobserved_md5\tmd5_status\toutput_file\n" \
    > "${VALIDATION_LOG}"

echo "Reading manifest: ${MANIFEST}"

# ------------------------------------------------------------------
# Download and checksum validation
# ------------------------------------------------------------------

tail -n +2 "${MANIFEST}" |
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
    [[ -z "${sample_id:-}" ]] && continue
    [[ -z "${file_url:-}" ]] && continue

    outfile="${OUTDIR}/${sample_id}.bam"
    partial_file="${outfile}.aria2"
    logfile="${LOGDIR}/${sample_id}.download.log"

    echo "----------------------------------------"
    echo "Sample: ${sample_id}"
    echo "Study: ${study_accession}"
    echo "Cancer type: ${cancer_type}"
    echo "Assay type: ${assay_type}"
    echo "Source: ${source}"
    echo "URL: ${file_url}"
    echo "Output: ${outfile}"
    echo "Log: ${logfile}"

    # Remove incomplete transfer state before restarting.
    if [[ -f "${partial_file}" ]]; then
        echo "Incomplete download state detected for ${sample_id}."
        echo "Removing partial BAM and aria2 control file."
        rm -f "${outfile}" "${partial_file}"
    fi

    # Download only when the BAM is not already present.
    if [[ -f "${outfile}" ]]; then
        echo "File already exists; proceeding to checksum validation."
    else
        echo "Downloading ${sample_id}..."

        if ! aria2c \
            -x 8 \
            -s 8 \
            -k 1M \
            -c \
            --file-allocation=none \
            -o "${sample_id}.bam" \
            -d "${OUTDIR}" \
            "${file_url}" |
            tee "${logfile}"
        then
            echo "ERROR: Download command failed for ${sample_id}."

            printf \
                "%s\t%s\t%s\t%s\tNA\tFAILED_DOWNLOAD\taria2c exited with an error\n" \
                "${sample_id}" \
                "${study_accession}" \
                "${file_url}" \
                "${expected_md5}" \
                >> "${FAIL_LOG}"

            continue
        fi
    fi

    if [[ ! -f "${outfile}" ]]; then
        echo "ERROR: Output BAM is missing for ${sample_id}."

        printf \
            "%s\t%s\t%s\t%s\tNA\tFAILED_DOWNLOAD\tOutput BAM missing after download\n" \
            "${sample_id}" \
            "${study_accession}" \
            "${file_url}" \
            "${expected_md5}" \
            >> "${FAIL_LOG}"

        continue
    fi

    if [[ -f "${partial_file}" ]]; then
        echo "ERROR: Download remains incomplete for ${sample_id}."

        printf \
            "%s\t%s\t%s\t%s\tNA\tFAILED_INCOMPLETE\taria2 control file remains after download\n" \
            "${sample_id}" \
            "${study_accession}" \
            "${file_url}" \
            "${expected_md5}" \
            >> "${FAIL_LOG}"

        continue
    fi

    echo "Calculating MD5 checksum for ${sample_id}..."
    observed_md5="$(md5sum "${outfile}" | awk '{print $1}')"

    if [[ "${observed_md5}" == "${expected_md5}" ]]; then
        echo "MD5 PASS: ${sample_id}"

        printf \
            "%s\t%s\t%s\t%s\t%s\tPASS\t%s\n" \
            "${sample_id}" \
            "${study_accession}" \
            "${file_url}" \
            "${expected_md5}" \
            "${observed_md5}" \
            "${outfile}" \
            >> "${VALIDATION_LOG}"
    else
        echo "MD5 FAIL: ${sample_id}"

        quarantine_file="${QUARANTINE_DIR}/${sample_id}.bam"

        # Avoid silently replacing an existing quarantined file.
        if [[ -e "${quarantine_file}" ]]; then
            timestamp="$(date +%Y%m%d_%H%M%S)"
            quarantine_file="${QUARANTINE_DIR}/${sample_id}_${timestamp}.bam"
        fi

        mv "${outfile}" "${quarantine_file}"

        printf \
            "%s\t%s\t%s\t%s\t%s\tFAILED_MD5\tMoved to %s\n" \
            "${sample_id}" \
            "${study_accession}" \
            "${file_url}" \
            "${expected_md5}" \
            "${observed_md5}" \
            "${quarantine_file}" \
            >> "${FAIL_LOG}"

        printf \
            "%s\t%s\t%s\t%s\t%s\tFAIL\t%s\n" \
            "${sample_id}" \
            "${study_accession}" \
            "${file_url}" \
            "${expected_md5}" \
            "${observed_md5}" \
            "${quarantine_file}" \
            >> "${VALIDATION_LOG}"
    fi
done

echo "Download and checksum-validation step complete."
echo "Validated source-BAM directory: ${OUTDIR}"
echo "Validation report: ${VALIDATION_LOG}"
echo "Failure report: ${FAIL_LOG}"
echo "Quarantine directory: ${QUARANTINE_DIR}"