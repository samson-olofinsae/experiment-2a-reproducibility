#!/usr/bin/env python3

"""
Generate benchmarking tables for the baseline concordance experiment.

This script parses the raw outputs produced by:

    scripts/compute_flagstat.sh
    scripts/compute_idxstats.sh

It generates:

    results/metrics_table.tsv
    results/discrepancy_summary.tsv

The analysis compares corresponding mapped-read quantities reported by
samtools flagstat and samtools idxstats across the deterministic nested
downsampling conditions used in the baseline concordance experiment.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError as exc:
    raise SystemExit(
        "ERROR: pandas is required but is not installed.\n"
        "Create the project environment from environment.yml or install pandas."
    ) from exc


FLAGSTAT_DIR = Path("results/flagstat")
IDXSTATS_DIR = Path("results/idxstats")
DOWNSAMPLED_DIR = Path("data/downsampled")
RESULTS_DIR = Path("results")

METRICS_TABLE = RESULTS_DIR / "metrics_table.tsv"
SUMMARY_TABLE = RESULTS_DIR / "discrepancy_summary.tsv"

EXPECTED_SAMPLES = {
    "SRR10353673",
    "SRR10353687",
    "SRR10353710",
}

EXPECTED_DEPTHS = {
    10,
    20,
    30,
    40,
    50,
    60,
    70,
    80,
    90,
}

EXPECTED_COMPARISONS = len(EXPECTED_SAMPLES) * len(EXPECTED_DEPTHS)

FILE_PATTERN = re.compile(
    r"(?P<sample>.+?)_(?P<depth_percent>\d+)pct\."
    r"(?P<metric_type>flagstat|idxstats)\.txt$"
)


def parse_flagstat(path: Path) -> tuple[int, int]:
    """
    Parse total-read and mapped-read counts from raw samtools flagstat output.

    Both QC-passed and QC-failed counts are summed so that F represents the
    complete mapped category used in the baseline concordance analysis.
    """
    total: int | None = None
    mapped: int | None = None

    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if "in total" in line:
                    parts = line.split()

                    if len(parts) < 3:
                        raise ValueError(
                            f"Malformed total-read line in flagstat file: {path}"
                        )

                    total = int(parts[0]) + int(parts[2])

                elif " mapped (" in line and mapped is None:
                    parts = line.split()

                    if len(parts) < 3:
                        raise ValueError(
                            f"Malformed mapped-read line in flagstat file: {path}"
                        )

                    mapped = int(parts[0]) + int(parts[2])

    except OSError as exc:
        raise OSError(f"Could not read flagstat file: {path}") from exc

    if total is None or mapped is None:
        raise ValueError(
            f"Could not parse total and mapped counts from flagstat file: {path}"
        )

    return total, mapped


def parse_idxstats(path: Path) -> tuple[int, int]:
    """
    Parse total-read and mapped-read counts from raw samtools idxstats output.

    I is defined as the sum of mapped alignment counts reported across all
    reference-sequence rows. Total reads are calculated as the sum of mapped
    and unmapped counts across the parsed idxstats rows.
    """
    total_mapped = 0
    total_unmapped = 0
    parsed_rows = 0

    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                parts = line.rstrip("\n").split("\t")

                if len(parts) < 4:
                    if line.strip():
                        print(
                            "WARNING: Skipping malformed idxstats line "
                            f"{line_number} in {path}",
                            file=sys.stderr,
                        )
                    continue

                try:
                    mapped = int(parts[2])
                    unmapped = int(parts[3])
                except ValueError as exc:
                    raise ValueError(
                        f"Non-integer mapped or unmapped count on line "
                        f"{line_number} of {path}"
                    ) from exc

                total_mapped += mapped
                total_unmapped += unmapped
                parsed_rows += 1

    except OSError as exc:
        raise OSError(f"Could not read idxstats file: {path}") from exc

    if parsed_rows == 0:
        raise ValueError(f"No valid idxstats rows were parsed from: {path}")

    return total_mapped + total_unmapped, total_mapped


def validate_input_directories() -> None:
    """Confirm that the expected raw-metric directories are available."""
    missing_directories = [
        str(path)
        for path in (FLAGSTAT_DIR, IDXSTATS_DIR)
        if not path.is_dir()
    ]

    if missing_directories:
        missing_text = "\n".join(
            f"  - {path}" for path in missing_directories
        )

        raise SystemExit(
            "ERROR: Required input directories were not found:\n"
            f"{missing_text}\n"
            "Run the metric-extraction workflow before generating the tables."
        )


def collect_rows() -> list[dict[str, object]]:
    """Parse all paired flagstat and idxstats outputs."""
    flagstat_files = sorted(FLAGSTAT_DIR.glob("*.flagstat.txt"))

    if not flagstat_files:
        raise SystemExit(
            f"ERROR: No flagstat files were found in {FLAGSTAT_DIR}."
        )

    rows: list[dict[str, object]] = []
    observed_keys: set[tuple[str, int]] = set()

    for flagstat_path in flagstat_files:
        match = FILE_PATTERN.fullmatch(flagstat_path.name)

        if not match:
            raise ValueError(
                "Unexpected flagstat filename. Expected the pattern "
                "'<sample>_<depth>pct.flagstat.txt': "
                f"{flagstat_path.name}"
            )

        sample = match.group("sample")
        depth_percent = int(match.group("depth_percent"))
        depth_fraction = depth_percent / 100

        key = (sample, depth_percent)

        if key in observed_keys:
            raise ValueError(
                "Duplicate sample-depth combination encountered: "
                f"{sample}, {depth_percent}%"
            )

        observed_keys.add(key)

        idxstats_path = (
            IDXSTATS_DIR
            / f"{sample}_{depth_percent}pct.idxstats.txt"
        )

        downsampled_bam = (
            DOWNSAMPLED_DIR
            / f"{sample}_{depth_percent}pct.bam"
        )

        if not idxstats_path.is_file():
            raise FileNotFoundError(
                f"Missing matching idxstats file: {idxstats_path}"
            )

        flagstat_total, flagstat_mapped = parse_flagstat(flagstat_path)
        idxstats_total, idxstats_mapped = parse_idxstats(idxstats_path)

        signed_discrepancy = flagstat_mapped - idxstats_mapped
        absolute_discrepancy = abs(signed_discrepancy)

        percentage_discrepancy = (
            absolute_discrepancy / flagstat_mapped * 100
            if flagstat_mapped > 0
            else None
        )

        exact_concordance = (
            flagstat_mapped == idxstats_mapped
        )

        rows.append(
            {
                "sample": sample,
                "depth_fraction": depth_fraction,
                "depth_percent": depth_percent,
                "downsampled_bam": str(downsampled_bam),
                "flagstat_total": flagstat_total,
                "flagstat_mapped": flagstat_mapped,
                "idxstats_total": idxstats_total,
                "idxstats_mapped": idxstats_mapped,
                "signed_discrepancy": signed_discrepancy,
                "absolute_discrepancy": absolute_discrepancy,
                "percentage_discrepancy": (
                    round(percentage_discrepancy, 6)
                    if percentage_discrepancy is not None
                    else None
                ),
                "exact_concordance": exact_concordance,
            }
        )

    return rows


def validate_experiment_design(
    metrics: pd.DataFrame,
) -> None:
    """
    Confirm that the parsed outputs represent the complete baseline design.
    """
    if len(metrics) != EXPECTED_COMPARISONS:
        raise ValueError(
            "Unexpected number of benchmarking comparisons: "
            f"found {len(metrics)}, expected {EXPECTED_COMPARISONS}."
        )

    observed_samples = set(
        metrics["sample"].astype(str)
    )

    if observed_samples != EXPECTED_SAMPLES:
        raise ValueError(
            "Unexpected baseline source-BAM set.\n"
            f"Observed: {sorted(observed_samples)}\n"
            f"Expected: {sorted(EXPECTED_SAMPLES)}"
        )

    duplicate_mask = metrics.duplicated(
        subset=[
            "sample",
            "depth_percent",
        ],
        keep=False,
    )

    if duplicate_mask.any():
        duplicate_rows = metrics.loc[
            duplicate_mask,
            [
                "sample",
                "depth_percent",
            ],
        ]

        raise ValueError(
            "Duplicate sample-depth combinations were found:\n"
            f"{duplicate_rows.to_string(index=False)}"
        )

    for sample in sorted(EXPECTED_SAMPLES):
        sample_depths = set(
            metrics.loc[
                metrics["sample"] == sample,
                "depth_percent",
            ].astype(int)
        )

        if sample_depths != EXPECTED_DEPTHS:
            missing_depths = sorted(
                EXPECTED_DEPTHS - sample_depths
            )

            unexpected_depths = sorted(
                sample_depths - EXPECTED_DEPTHS
            )

            raise ValueError(
                f"Incomplete depth series for sample {sample}.\n"
                f"Missing depths: {missing_depths}\n"
                f"Unexpected depths: {unexpected_depths}"
            )


def summarise_group(
    group: pd.DataFrame,
    *,
    summary_level: str,
    depth_percent: int | None,
) -> dict[str, object]:
    """
    Summarise concordance and discrepancy outcomes for a set of observations.
    """
    n_observations = len(group)

    if n_observations == 0:
        raise ValueError(
            "Cannot summarise an empty observation group."
        )

    exact_count = int(
        group["exact_concordance"].sum()
    )

    percentage_values = (
        group["percentage_discrepancy"]
        .dropna()
    )

    max_percentage_discrepancy = (
        float(percentage_values.max())
        if not percentage_values.empty
        else None
    )

    return {
        "summary_level": summary_level,
        "depth_percent": depth_percent,
        "n_observations": n_observations,
        "mean_flagstat_mapped": round(
            float(group["flagstat_mapped"].mean()),
            2,
        ),
        "mean_idxstats_mapped": round(
            float(group["idxstats_mapped"].mean()),
            2,
        ),
        "exact_concordance_count": exact_count,
        "exact_concordance_pct": round(
            exact_count / n_observations * 100,
            6,
        ),
        "max_absolute_discrepancy": int(
            group["absolute_discrepancy"].max()
        ),
        "max_percentage_discrepancy": (
            round(
                max_percentage_discrepancy,
                6,
            )
            if max_percentage_discrepancy is not None
            else None
        ),
    }


def build_summary(
    metrics: pd.DataFrame,
) -> pd.DataFrame:
    """
    Generate depth-level and overall baseline concordance summaries.
    """
    summary_rows: list[dict[str, object]] = []

    for depth_percent in sorted(EXPECTED_DEPTHS):
        depth_group = metrics.loc[
            metrics["depth_percent"] == depth_percent
        ]

        summary_rows.append(
            summarise_group(
                depth_group,
                summary_level="depth",
                depth_percent=depth_percent,
            )
        )

    summary_rows.append(
        summarise_group(
            metrics,
            summary_level="overall",
            depth_percent=None,
        )
    )

    summary = pd.DataFrame(summary_rows)

    summary["depth_percent"] = (
        summary["depth_percent"].astype("Int64")
    )

    return summary


def validate_results(
    metrics: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    """
    Validate internally derived concordance and discrepancy quantities.
    """
    expected_absolute = (
        metrics["signed_discrepancy"].abs()
    )

    if not metrics[
        "absolute_discrepancy"
    ].equals(expected_absolute):
        raise ValueError(
            "Absolute discrepancy does not match "
            "the absolute signed discrepancy."
        )

    expected_exact = (
        metrics["flagstat_mapped"]
        == metrics["idxstats_mapped"]
    )

    if not metrics[
        "exact_concordance"
    ].equals(expected_exact):
        raise ValueError(
            "Exact-concordance indicators are inconsistent "
            "with mapped-read counts."
        )

    if (
        metrics["absolute_discrepancy"] < 0
    ).any():
        raise ValueError(
            "Absolute discrepancies cannot be negative."
        )

    non_null_percentage = (
        metrics["percentage_discrepancy"]
        .dropna()
    )

    if (
        non_null_percentage < 0
    ).any():
        raise ValueError(
            "Percentage discrepancies cannot be negative."
        )

    nonpositive_flagstat = (
        metrics["flagstat_mapped"] <= 0
    )

    if metrics.loc[
        nonpositive_flagstat,
        "percentage_discrepancy",
    ].notna().any():
        raise ValueError(
            "Percentage discrepancy must be undefined "
            "when flagstat_mapped is not positive."
        )

    overall_rows = summary.loc[
        summary["summary_level"] == "overall"
    ]

    if len(overall_rows) != 1:
        raise ValueError(
            "Expected exactly one overall summary row."
        )

    overall = overall_rows.iloc[0]

    if (
        int(overall["n_observations"])
        != EXPECTED_COMPARISONS
    ):
        raise ValueError(
            "Overall summary does not contain the expected "
            "number of paired observations."
        )

    if (
        int(overall["exact_concordance_count"])
        != int(metrics["exact_concordance"].sum())
    ):
        raise ValueError(
            "Overall exact-concordance count is inconsistent "
            "with the observation-level table."
        )

    if (
        int(overall["max_absolute_discrepancy"])
        != int(metrics["absolute_discrepancy"].max())
    ):
        raise ValueError(
            "Overall maximum absolute discrepancy is inconsistent "
            "with the observation-level table."
        )


def main() -> None:
    """
    Generate and validate the canonical baseline concordance tables.
    """
    validate_input_directories()

    RESULTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    rows = collect_rows()

    if not rows:
        raise SystemExit(
            "ERROR: No valid paired metric records were generated."
        )

    metrics = (
        pd.DataFrame(rows)
        .sort_values(
            [
                "sample",
                "depth_percent",
            ],
            kind="stable",
        )
        .reset_index(drop=True)
    )

    validate_experiment_design(metrics)

    summary = build_summary(metrics)

    validate_results(
        metrics,
        summary,
    )

    metrics.to_csv(
        METRICS_TABLE,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary.to_csv(
        SUMMARY_TABLE,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    overall = summary.loc[
        summary["summary_level"] == "overall"
    ].iloc[0]

    maximum_percentage = (
        overall["max_percentage_discrepancy"]
    )

    if pd.isna(maximum_percentage):
        percentage_text = "NA"
    else:
        percentage_text = (
            f"{float(maximum_percentage):.6f}%"
        )

    print(
        f"Wrote {METRICS_TABLE} with {len(metrics)} paired observations "
        f"({len(EXPECTED_SAMPLES)} source BAMs × "
        f"{len(EXPECTED_DEPTHS)} nested depth conditions)."
    )

    print(
        f"Wrote {SUMMARY_TABLE} with "
        f"{len(summary)} summary rows."
    )

    print(
        "Baseline concordance summary: "
        f"{int(overall['exact_concordance_count'])}/"
        f"{int(overall['n_observations'])} exact; "
        "maximum absolute discrepancy = "
        f"{int(overall['max_absolute_discrepancy'])}; "
        "maximum percentage discrepancy = "
        f"{percentage_text}."
    )

    print(
        "Baseline concordance table generation "
        "completed successfully."
    )


if __name__ == "__main__":
    try:
        main()
    except (
        FileNotFoundError,
        OSError,
        ValueError,
    ) as exc:
        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc