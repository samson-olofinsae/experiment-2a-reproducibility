#!/usr/bin/env python3

"""
Generate publication figures and figure-level summary statistics for the
baseline mapped-read concordance experiment.

Input:
    results/metrics_table.tsv

Outputs:
    results/figures/figure2_mean_mapped_reads_by_depth.png
    results/figures/figure2_mean_mapped_reads_by_depth.pdf

    results/figures/figure3_signed_discrepancy_by_depth.png
    results/figures/figure3_signed_discrepancy_by_depth.pdf

    results/figures/figure_summary_stats.txt

Figure 2 summarises mapped-read counts reported by samtools flagstat
and samtools idxstats across deterministic nested downsampling conditions.

Figure 3 displays all 27 baseline observations individually as an
observation matrix, avoiding graphical overlap when signed discrepancies
are identical.

Pearson correlation is intentionally not used because correlation does not
measure equality or agreement between paired mapped-read quantities.
"""

from __future__ import annotations

import sys
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.ticker import FuncFormatter
except ImportError as exc:
    raise SystemExit(
        "ERROR: pandas and matplotlib are required but are not installed.\n"
        "Create the project environment from environment.yml."
    ) from exc


RESULTS_DIR = Path("results")
METRICS_TABLE = RESULTS_DIR / "metrics_table.tsv"
FIGURES_DIR = RESULTS_DIR / "figures"

EXPECTED_ROWS = 27

EXPECTED_SAMPLES = [
    "SRR10353673",
    "SRR10353687",
    "SRR10353710",
]

EXPECTED_DEPTHS = [
    10,
    20,
    30,
    40,
    50,
    60,
    70,
    80,
    90,
]

REQUIRED_COLUMNS = {
    "sample",
    "depth_fraction",
    "depth_percent",
    "downsampled_bam",
    "flagstat_total",
    "flagstat_mapped",
    "idxstats_total",
    "idxstats_mapped",
    "signed_discrepancy",
    "absolute_discrepancy",
    "percentage_discrepancy",
    "exact_concordance",
}

FIGURE_FILES = {
    "figure2_png": (
        FIGURES_DIR / "figure2_mean_mapped_reads_by_depth.png"
    ),
    "figure2_pdf": (
        FIGURES_DIR / "figure2_mean_mapped_reads_by_depth.pdf"
    ),
    "figure3_png": (
        FIGURES_DIR / "figure3_signed_discrepancy_by_depth.png"
    ),
    "figure3_pdf": (
        FIGURES_DIR / "figure3_signed_discrepancy_by_depth.pdf"
    ),
    "summary": (
        FIGURES_DIR / "figure_summary_stats.txt"
    ),
}

LEGACY_FIGURE_FILES = [
    FIGURES_DIR / "figure3_flagstat_vs_idxstats_agreement.png",
    FIGURES_DIR / "figure3_flagstat_vs_idxstats_agreement.pdf",
    FIGURES_DIR / "figure4_mapped_difference_pct_by_depth.png",
    FIGURES_DIR / "figure4_mapped_difference_pct_by_depth.pdf",
    FIGURES_DIR / "figure5_sample_mapped_reads_by_depth.png",
    FIGURES_DIR / "figure5_sample_mapped_reads_by_depth.pdf",
]


def reads_in_millions(
    value: float,
    position: float,
) -> str:
    """Format mapped-read counts in millions on plot axes."""
    del position
    return f"{value / 1_000_000:.0f}"


def load_and_validate_metrics() -> pd.DataFrame:
    """Load and validate the canonical baseline metrics table."""
    if not METRICS_TABLE.is_file():
        raise FileNotFoundError(
            f"Metrics table not found: {METRICS_TABLE}\n"
            "Run scripts/generate_metrics_table.py first."
        )

    try:
        metrics = pd.read_csv(
            METRICS_TABLE,
            sep="\t",
        )
    except Exception as exc:
        raise ValueError(
            f"Could not read metrics table: {METRICS_TABLE}"
        ) from exc

    missing_columns = REQUIRED_COLUMNS - set(metrics.columns)

    if missing_columns:
        raise ValueError(
            "Metrics table is missing required columns: "
            f"{sorted(missing_columns)}"
        )

    if len(metrics) != EXPECTED_ROWS:
        raise ValueError(
            f"Unexpected number of rows in {METRICS_TABLE}: "
            f"found {len(metrics)}, expected {EXPECTED_ROWS}."
        )

    observed_samples = set(
        metrics["sample"].astype(str)
    )

    if observed_samples != set(EXPECTED_SAMPLES):
        raise ValueError(
            "Unexpected source-BAM set in metrics table.\n"
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
            "Duplicate sample-depth combinations found:\n"
            f"{duplicate_rows.to_string(index=False)}"
        )

    for sample in EXPECTED_SAMPLES:
        observed_depths = set(
            metrics.loc[
                metrics["sample"] == sample,
                "depth_percent",
            ].astype(int)
        )

        if observed_depths != set(EXPECTED_DEPTHS):
            raise ValueError(
                f"Unexpected depth series for sample {sample}.\n"
                f"Observed: {sorted(observed_depths)}\n"
                f"Expected: {EXPECTED_DEPTHS}"
            )

    required_numeric_columns = [
        "depth_fraction",
        "depth_percent",
        "flagstat_total",
        "flagstat_mapped",
        "idxstats_total",
        "idxstats_mapped",
        "signed_discrepancy",
        "absolute_discrepancy",
    ]

    for column in required_numeric_columns:
        if metrics[column].isna().any():
            raise ValueError(
                f"Column contains missing values: {column}"
            )

    expected_signed = (
        metrics["flagstat_mapped"]
        - metrics["idxstats_mapped"]
    )

    if not metrics[
        "signed_discrepancy"
    ].equals(expected_signed):
        raise ValueError(
            "signed_discrepancy is inconsistent with "
            "flagstat_mapped - idxstats_mapped."
        )

    expected_absolute = (
        metrics["signed_discrepancy"].abs()
    )

    if not metrics[
        "absolute_discrepancy"
    ].equals(expected_absolute):
        raise ValueError(
            "absolute_discrepancy is inconsistent with "
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
            "exact_concordance is inconsistent with "
            "the mapped-read counts."
        )

    nonpositive_flagstat = (
        metrics["flagstat_mapped"] <= 0
    )

    if metrics.loc[
        nonpositive_flagstat,
        "percentage_discrepancy",
    ].notna().any():
        raise ValueError(
            "percentage_discrepancy must be undefined "
            "when flagstat_mapped is not positive."
        )

    return (
        metrics.sort_values(
            [
                "sample",
                "depth_percent",
            ],
            kind="stable",
        )
        .reset_index(drop=True)
    )


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


def remove_legacy_outputs() -> None:
    """Remove obsolete outputs from the retired plotting workflow."""
    for path in LEGACY_FIGURE_FILES:
        if path.exists():
            path.unlink()


def save_figure(
    png_path: Path,
    pdf_path: Path,
) -> None:
    """Save the current figure in PNG and PDF formats."""
    plt.tight_layout()

    plt.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    plt.close()

    for path in (
        png_path,
        pdf_path,
    ):
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise OSError(
                "Figure output was not created "
                f"or is empty: {path}"
            )


def build_depth_summary(
    metrics: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate mean mapped-read counts at each nested depth condition."""
    return (
        metrics.groupby(
            "depth_percent",
            sort=True,
        )[
            [
                "flagstat_mapped",
                "idxstats_mapped",
            ]
        ]
        .mean()
        .reset_index()
    )


def generate_figure2(
    depth_summary: pd.DataFrame,
) -> None:
    """Generate Figure 2: mean mapped-read counts across depth."""
    plt.figure(
        figsize=(7, 5)
    )

    plt.plot(
        depth_summary["depth_percent"],
        depth_summary["flagstat_mapped"],
        marker="o",
        linewidth=2.2,
        label="samtools flagstat",
    )

    plt.plot(
        depth_summary["depth_percent"],
        depth_summary["idxstats_mapped"],
        marker="s",
        linestyle="--",
        linewidth=1.9,
        label="samtools idxstats",
    )

    plt.xlabel(
        "Downsampling fraction (%)"
    )

    plt.ylabel(
        "Mean mapped-read count (millions)"
    )

    plt.gca().yaxis.set_major_formatter(
        FuncFormatter(
            reads_in_millions
        )
    )

    plt.grid(
        True,
        alpha=0.25,
    )

    plt.legend(
        frameon=True,
        loc="upper left",
    )

    save_figure(
        FIGURE_FILES["figure2_png"],
        FIGURE_FILES["figure2_pdf"],
    )


def generate_figure3(
    metrics: pd.DataFrame,
) -> None:
    """
    Generate Figure 3 as a sample-by-depth observation matrix.

    Every sample-depth combination is represented separately. Each marker
    contains the observed signed discrepancy:

        Δ = flagstat_mapped - idxstats_mapped

    A single marker appearance is used throughout because colour does not
    encode sample, depth, or discrepancy category in this figure.
    """
    fig, ax = plt.subplots(
        figsize=(8, 4.2)
    )

    y_positions = {
        sample: index
        for index, sample in enumerate(EXPECTED_SAMPLES)
    }

    for row in metrics.itertuples(index=False):
        x = int(row.depth_percent)
        y = y_positions[str(row.sample)]
        discrepancy = int(row.signed_discrepancy)

        ax.scatter(
            x,
            y,
            s=430,
            marker="o",
            facecolor="white",
            edgecolor="black",
            linewidth=1.4,
            zorder=2,
        )

        ax.text(
            x,
            y,
            str(discrepancy),
            ha="center",
            va="center",
            fontsize=9,
            color="black",
            zorder=3,
        )

    ax.set_xlabel(
        "Downsampling fraction (%)"
    )

    ax.set_ylabel(
        "Source BAM"
    )

    ax.set_xticks(
        EXPECTED_DEPTHS
    )

    ax.set_yticks(
        list(range(len(EXPECTED_SAMPLES)))
    )

    ax.set_yticklabels(
        EXPECTED_SAMPLES
    )

    ax.set_xlim(
        min(EXPECTED_DEPTHS) - 4,
        max(EXPECTED_DEPTHS) + 4,
    )

    ax.set_ylim(
        -0.6,
        len(EXPECTED_SAMPLES) - 0.4,
    )

    ax.grid(
        True,
        axis="x",
        alpha=0.18,
        zorder=0,
    )

    ax.text(
        0.5,
        -0.22,
        "Marker labels show signed discrepancy Δ = F − I (reads)",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10,
    )

    save_figure(
        FIGURE_FILES["figure3_png"],
        FIGURE_FILES["figure3_pdf"],
    )


def write_summary_statistics(
    metrics: pd.DataFrame,
) -> None:
    """Write figure-level baseline concordance statistics."""
    summary_path = FIGURE_FILES["summary"]

    comparison_count = len(metrics)

    exact_count = int(
        metrics["exact_concordance"].sum()
    )

    exact_percentage = (
        exact_count
        / comparison_count
        * 100
    )

    maximum_absolute = int(
        metrics["absolute_discrepancy"].max()
    )

    percentage_values = (
        metrics["percentage_discrepancy"]
        .dropna()
    )

    if percentage_values.empty:
        maximum_percentage = None
    else:
        maximum_percentage = float(
            percentage_values.max()
        )

    with summary_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "Baseline mapped-read concordance summary\n"
        )

        handle.write(
            f"Number of paired observations: "
            f"{comparison_count}\n"
        )

        handle.write(
            f"Exact concordance count: "
            f"{exact_count}\n"
        )

        handle.write(
            "Exact concordance percentage: "
            f"{exact_percentage:.6f}\n"
        )

        handle.write(
            "Maximum absolute mapped-read discrepancy: "
            f"{maximum_absolute}\n"
        )

        if maximum_percentage is None:
            handle.write(
                "Maximum percentage mapped-read discrepancy: "
                "NA\n"
            )
        else:
            handle.write(
                "Maximum percentage mapped-read discrepancy: "
                f"{maximum_percentage:.6f}\n"
            )

    if (
        not summary_path.is_file()
        or summary_path.stat().st_size == 0
    ):
        raise OSError(
            "Summary statistics output was not created: "
            f"{summary_path}"
        )


def main() -> None:
    """Generate and verify revised baseline figures."""
    metrics = load_and_validate_metrics()

    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    configure_plotting()

    remove_legacy_outputs()

    depth_summary = build_depth_summary(
        metrics
    )

    generate_figure2(
        depth_summary
    )

    generate_figure3(
        metrics
    )

    write_summary_statistics(
        metrics
    )

    exact_count = int(
        metrics["exact_concordance"].sum()
    )

    comparison_count = len(metrics)

    maximum_absolute = int(
        metrics["absolute_discrepancy"].max()
    )

    percentage_values = (
        metrics["percentage_discrepancy"]
        .dropna()
    )

    if percentage_values.empty:
        maximum_percentage_text = "NA"
    else:
        maximum_percentage_text = (
            f"{float(percentage_values.max()):.6f}%"
        )

    print(
        "Generated revised baseline publication figures "
        f"in {FIGURES_DIR}."
    )

    print(
        "Baseline concordance summary: "
        f"{exact_count}/{comparison_count} exact; "
        "maximum absolute discrepancy = "
        f"{maximum_absolute}; "
        "maximum percentage discrepancy = "
        f"{maximum_percentage_text}."
    )

    for path in FIGURE_FILES.values():
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