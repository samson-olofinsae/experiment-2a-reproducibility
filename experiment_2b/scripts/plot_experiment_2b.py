#!/usr/bin/env python3

"""
Generate the publication figure and summary statistics for the
stale-index positive-control experiment.

Input:
    experiment_2b/results/experiment_2b_results.tsv

Outputs:
    experiment_2b/results/figures/
        figure4_stale_index_expected_vs_observed_discrepancy.png
        figure4_stale_index_expected_vs_observed_discrepancy.pdf
        figure4_summary_stats.txt

Figure 4 displays expected and observed signed discrepancies separately
within each condition so that exact quantitative recovery remains visible
rather than being obscured by graphical overlap.
"""

from __future__ import annotations

import sys
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import pandas as pd
except ImportError as exc:
    raise SystemExit(
        "ERROR: pandas and matplotlib are required but are not installed.\n"
        "Create the project environment from environment.yml."
    ) from exc


EXPERIMENT_DIR = Path("experiment_2b")

RESULTS_TABLE = (
    EXPERIMENT_DIR
    / "results"
    / "experiment_2b_results.tsv"
)

FIGURES_DIR = (
    EXPERIMENT_DIR
    / "results"
    / "figures"
)

FIGURE_PNG = (
    FIGURES_DIR
    / "figure4_stale_index_expected_vs_observed_discrepancy.png"
)

FIGURE_PDF = (
    FIGURES_DIR
    / "figure4_stale_index_expected_vs_observed_discrepancy.pdf"
)

SUMMARY_FILE = (
    FIGURES_DIR
    / "figure4_summary_stats.txt"
)

EXPECTED_CONDITIONS = [
    "baseline",
    "PC1_k1000",
    "PC2_k10000",
]

DISPLAY_LABELS = {
    "baseline": "Baseline",
    "PC1_k1000": "PC1\nk = 1,000",
    "PC2_k10000": "PC2\nk = 10,000",
}

REQUIRED_COLUMNS = {
    "condition",
    "baseline_mapped_count",
    "removed_mapped_records",
    "flagstat_mapped",
    "idxstats_mapped",
    "delta_expected",
    "delta_observed",
    "residual",
}


def load_and_validate_results() -> pd.DataFrame:
    """Load and validate canonical positive-control results."""
    if not RESULTS_TABLE.is_file():
        raise FileNotFoundError(
            f"Results table not found: {RESULTS_TABLE}\n"
            "Run experiment_2b/scripts/run_experiment_2b.sh first."
        )

    try:
        results = pd.read_csv(
            RESULTS_TABLE,
            sep="\t",
        )
    except Exception as exc:
        raise ValueError(
            f"Could not read positive-control results: {RESULTS_TABLE}"
        ) from exc

    missing_columns = REQUIRED_COLUMNS - set(results.columns)

    if missing_columns:
        raise ValueError(
            "Positive-control results table is missing required columns: "
            f"{sorted(missing_columns)}"
        )

    if len(results) != len(EXPECTED_CONDITIONS):
        raise ValueError(
            "Unexpected number of positive-control rows: "
            f"found {len(results)}, "
            f"expected {len(EXPECTED_CONDITIONS)}."
        )

    observed_conditions = (
        results["condition"]
        .astype(str)
        .tolist()
    )

    if observed_conditions != EXPECTED_CONDITIONS:
        raise ValueError(
            "Unexpected positive-control condition order.\n"
            f"Observed: {observed_conditions}\n"
            f"Expected: {EXPECTED_CONDITIONS}"
        )

    numeric_columns = [
        "baseline_mapped_count",
        "removed_mapped_records",
        "flagstat_mapped",
        "idxstats_mapped",
        "delta_expected",
        "delta_observed",
        "residual",
    ]

    for column in numeric_columns:
        if results[column].isna().any():
            raise ValueError(
                f"Column contains missing values: {column}"
            )

    calculated_observed = (
        results["flagstat_mapped"]
        - results["idxstats_mapped"]
    )

    if not results[
        "delta_observed"
    ].equals(calculated_observed):
        raise ValueError(
            "delta_observed is inconsistent with "
            "flagstat_mapped - idxstats_mapped."
        )

    calculated_residual = (
        results["delta_observed"]
        - results["delta_expected"]
    )

    if not results[
        "residual"
    ].equals(calculated_residual):
        raise ValueError(
            "residual is inconsistent with "
            "delta_observed - delta_expected."
        )

    baseline = results.iloc[0]

    if (
        baseline["removed_mapped_records"] != 0
        or baseline["delta_expected"] != 0
    ):
        raise ValueError(
            "Baseline row does not represent the expected unperturbed state."
        )

    for row in results.iloc[1:].itertuples(index=False):
        expected_delta = -int(
            row.removed_mapped_records
        )

        if int(row.delta_expected) != expected_delta:
            raise ValueError(
                "Prespecified discrepancy is inconsistent with "
                f"removed mapped records for condition {row.condition}."
            )

    return results


def configure_plotting() -> None:
    """Apply a consistent publication-style plotting configuration."""
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "figure.dpi": 300,
        }
    )


def save_figure() -> None:
    """Save Figure 4 in PNG and PDF formats."""
    plt.tight_layout()

    plt.savefig(
        FIGURE_PNG,
        dpi=300,
        bbox_inches="tight",
    )

    plt.savefig(
        FIGURE_PDF,
        bbox_inches="tight",
    )

    plt.close()

    for path in (
        FIGURE_PNG,
        FIGURE_PDF,
    ):
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise OSError(
                "Figure output was not created "
                f"or is empty: {path}"
            )


def generate_figure(
    results: pd.DataFrame,
) -> None:
    """
    Generate Figure 4 using horizontally separated expected and observed
    markers within each condition.

    Horizontal separation preserves visibility when expected and observed
    discrepancies are numerically identical. Neutral connector lines
    indicate pairing without encoding an additional variable.
    """
    base_positions = list(
        range(len(results))
    )

    offset = 0.12

    expected_positions = [
        position - offset
        for position in base_positions
    ]

    observed_positions = [
        position + offset
        for position in base_positions
    ]

    expected = (
        results["delta_expected"]
        .astype(int)
        .tolist()
    )

    observed = (
        results["delta_observed"]
        .astype(int)
        .tolist()
    )

    labels = [
        DISPLAY_LABELS[condition]
        for condition in results["condition"]
    ]

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    for x_expected, x_observed, y_expected, y_observed in zip(
        expected_positions,
        observed_positions,
        expected,
        observed,
    ):
        ax.plot(
            [
                x_expected,
                x_observed,
            ],
            [
                y_expected,
                y_observed,
            ],
            color="0.55",
            linewidth=1.0,
            zorder=1,
        )

    ax.scatter(
        expected_positions,
        expected,
        marker="o",
        s=75,
        label="Expected Δ",
        zorder=3,
    )

    ax.scatter(
        observed_positions,
        observed,
        marker="s",
        s=75,
        label="Observed Δ",
        zorder=3,
    )

    ax.axhline(
        0,
        color="0.45",
        linestyle=":",
        linewidth=1.0,
        zorder=0,
    )

    ax.set_xticks(
        base_positions
    )

    ax.set_xticklabels(
        labels
    )

    ax.set_xlabel(
        "Positive-control condition"
    )

    ax.set_ylabel(
        "Signed mapped-read discrepancy, Δ (reads)"
    )

    ax.grid(
        True,
        axis="y",
        alpha=0.22,
    )

    ax.legend(
        frameon=True,
        loc="lower left",
    )

    for index, (
        x_position,
        value,
    ) in enumerate(
        zip(
            observed_positions,
            observed,
        )
    ):
        if index == 0:
            ax.annotate(
                "0",
                (
                    x_position,
                    value,
                ),
                xytext=(8, -8),
                textcoords="offset points",
                ha="left",
                va="top",
                fontsize=9,
            )
        else:
            ax.annotate(
                f"{value:,}",
                (
                    x_position,
                    value,
                ),
                xytext=(7, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=9,
            )

    save_figure()


def write_summary_statistics(
    results: pd.DataFrame,
) -> None:
    """Write concise positive-control recovery statistics."""
    residuals = (
        results["residual"]
        .astype(int)
    )

    exact_recovery_count = int(
        residuals.eq(0).sum()
    )

    maximum_absolute_residual = int(
        residuals.abs().max()
    )

    with SUMMARY_FILE.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "Stale-index positive-control summary\n"
        )

        handle.write(
            f"Number of evaluated conditions: "
            f"{len(results)}\n"
        )

        handle.write(
            "Exact discrepancy-recovery count: "
            f"{exact_recovery_count}\n"
        )

        handle.write(
            "Maximum absolute residual: "
            f"{maximum_absolute_residual}\n"
        )

        for row in results.itertuples(
            index=False
        ):
            handle.write(
                f"{row.condition}: "
                f"removed_mapped_records="
                f"{int(row.removed_mapped_records)}, "
                f"delta_expected="
                f"{int(row.delta_expected)}, "
                f"delta_observed="
                f"{int(row.delta_observed)}, "
                f"residual="
                f"{int(row.residual)}\n"
            )

    if (
        not SUMMARY_FILE.is_file()
        or SUMMARY_FILE.stat().st_size == 0
    ):
        raise OSError(
            "Summary output was not created: "
            f"{SUMMARY_FILE}"
        )


def main() -> None:
    """Generate and verify the positive-control publication figure."""
    results = load_and_validate_results()

    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    configure_plotting()

    generate_figure(
        results
    )

    write_summary_statistics(
        results
    )

    exact_recovery_count = int(
        results["residual"]
        .eq(0)
        .sum()
    )

    maximum_absolute_residual = int(
        results["residual"]
        .abs()
        .max()
    )

    print(
        "Generated stale-index positive-control "
        f"publication figure in {FIGURES_DIR}."
    )

    print(
        "Positive-control recovery summary: "
        f"{exact_recovery_count}/{len(results)} conditions "
        "recovered exactly; "
        "maximum absolute residual = "
        f"{maximum_absolute_residual}."
    )

    for path in (
        FIGURE_PNG,
        FIGURE_PDF,
        SUMMARY_FILE,
    ):
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise OSError(
                "Expected output was not created "
                f"or is empty: {path}"
            )

        print(
            f"Verified output: {path}"
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