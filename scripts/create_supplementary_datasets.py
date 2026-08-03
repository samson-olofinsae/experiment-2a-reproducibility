#!/usr/bin/env python3

"""
Create Supplementary Dataset Excel files for Experiment 2a.

Supplementary Dataset S1
    metrics_table.tsv

Supplementary Dataset S2
    discrepancy_summary.tsv

These files are direct Excel representations of the canonical TSV outputs
used in the manuscript.
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

EXPECTED = ROOT / "expected_results"
SUPPLEMENTARY = ROOT / "supplementary"

SUPPLEMENTARY.mkdir(exist_ok=True)

FILES = [
    (
        "metrics_table.tsv",
        "Supplementary_Dataset_S1.xlsx",
        "Metrics_Table",
    ),
    (
        "discrepancy_summary.tsv",
        "Supplementary_Dataset_S2.xlsx",
        "Summary_Statistics",
    ),
]


for tsv_name, excel_name, sheet in FILES:

    tsv = EXPECTED / tsv_name

    if not tsv.exists():
        raise FileNotFoundError(tsv)

    df = pd.read_csv(tsv, sep="\t")

    out = SUPPLEMENTARY / excel_name

    with pd.ExcelWriter(out, engine="openpyxl") as writer:

        df.to_excel(
            writer,
            sheet_name=sheet,
            index=False,
        )

        metadata = pd.DataFrame(
            {
                "Field": [
                    "Experiment",
                    "Repository",
                    "Purpose",
                    "Generated from",
                ],
                "Value": [
                    "Experiment 2a",
                    "experiment-2a-reproducibility",
                    "Supplementary dataset accompanying the manuscript",
                    tsv_name,
                ],
            }
        )

        metadata.to_excel(
            writer,
            sheet_name="Metadata",
            index=False,
        )

    print(f"Wrote {out}")

print("\nFinished.")