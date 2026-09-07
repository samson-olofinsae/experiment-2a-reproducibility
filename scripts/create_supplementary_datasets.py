#!/usr/bin/env python3

"""
Create the supplementary Excel datasets for the baseline concordance
experiment.

Supplementary Dataset S1
    expected_results/metrics_table.tsv

Supplementary Dataset S2
    expected_results/discrepancy_summary.tsv

The Excel workbooks are direct representations of the canonical TSV outputs
used to support the revised manuscript. No additional statistical analysis is
performed by this script.
"""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent

EXPECTED = ROOT / "expected_results"
SUPPLEMENTARY = ROOT / "supplementary"

SUPPLEMENTARY.mkdir(
    exist_ok=True
)

MANUSCRIPT_TITLE = (
    "A Reproducible Computational Benchmarking Framework for Evaluating "
    "Concordance Between BAM-Derived Mapped-Read Counts"
)

FILES = [
    (
        "metrics_table.tsv",
        "Supplementary_Dataset_S1.xlsx",
        "Metrics_Table",
        (
            "Observation-level mapped-read metrics and discrepancy measures "
            "for the 27 paired baseline observations."
        ),
    ),
    (
        "discrepancy_summary.tsv",
        "Supplementary_Dataset_S2.xlsx",
        "Summary_Statistics",
        (
            "Depth-level and overall summaries of mapped-read concordance "
            "for the baseline experiment."
        ),
    ),
]


def write_workbook(
    tsv_name: str,
    excel_name: str,
    sheet_name: str,
    purpose: str,
) -> None:
    """Create one supplementary Excel workbook from a canonical TSV file."""
    tsv_path = EXPECTED / tsv_name

    if not tsv_path.is_file():
        raise FileNotFoundError(
            f"Canonical TSV file not found: {tsv_path}"
        )

    dataframe = pd.read_csv(
        tsv_path,
        sep="\t",
    )

    if dataframe.empty:
        raise ValueError(
            f"Canonical TSV file contains no data rows: {tsv_path}"
        )

    output_path = (
        SUPPLEMENTARY / excel_name
    )

    metadata = pd.DataFrame(
        {
            "Field": [
                "Study component",
                "Manuscript title",
                "Repository",
                "Purpose",
                "Generated from",
            ],
            "Value": [
                "Baseline concordance experiment",
                MANUSCRIPT_TITLE,
                "experiment-2a-reproducibility",
                purpose,
                tsv_name,
            ],
        }
    )

    with pd.ExcelWriter(
        output_path,
        engine="openpyxl",
    ) as writer:
        dataframe.to_excel(
            writer,
            sheet_name=sheet_name,
            index=False,
        )

        metadata.to_excel(
            writer,
            sheet_name="Metadata",
            index=False,
        )

    if (
        not output_path.is_file()
        or output_path.stat().st_size == 0
    ):
        raise OSError(
            f"Supplementary workbook was not created: {output_path}"
        )

    print(
        f"Wrote {output_path}"
    )


def main() -> None:
    """Generate and verify all baseline supplementary datasets."""
    for (
        tsv_name,
        excel_name,
        sheet_name,
        purpose,
    ) in FILES:
        write_workbook(
            tsv_name=tsv_name,
            excel_name=excel_name,
            sheet_name=sheet_name,
            purpose=purpose,
        )


if __name__ == "__main__":
    main()