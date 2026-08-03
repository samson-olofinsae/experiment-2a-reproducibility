#!/usr/bin/env python3

"""
Generate the publication figures and figure-level summary statistics
for Experiment 2a.

Input:
    results/metrics_table.tsv

Outputs:
    results/figures/figure2_mean_mapped_reads_by_depth.png
    results/figures/figure2_mean_mapped_reads_by_depth.pdf

    results/figures/figure3_flagstat_vs_idxstats_agreement.png
    results/figures/figure3_flagstat_vs_idxstats_agreement.pdf

    results/figures/figure4_mapped_difference_pct_by_depth.png
    results/figures/figure4_mapped_difference_pct_by_depth.pdf

    results/figures/figure5_sample_mapped_reads_by_depth.png
    results/figures/figure5_sample_mapped_reads_by_depth.pdf

    results/figures/figure_summary_stats.txt
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

REQUIRED_COLUMNS = {
    "sample",
    "depth_fraction",
    "depth_percent",
    "downsampled_bam",
    "flagstat_total",
    "flagstat_mapped",
    "idxstats_total",
    "idxstats_mapped",
    "mapped_diff",
    "mapped_diff_pct",
}

FIGURE_FILES = {
    "figure2_png": (
        FIGURES_DIR / "figure2_mean_mapped_reads_by_depth.png"
    ),
    "figure2_pdf": (
        FIGURES_DIR / "figure2_mean_mapped_reads_by_depth.pdf"
    ),
    "figure3_png": (
        FIGURES_DIR / "figure3_flagstat_vs_idxstats_agreement.png"
    ),
    "figure3_pdf": (
        FIGURES_DIR / "figure3_flagstat_vs_idxstats_agreement.pdf"
    ),
    "figure4_png": (
        FIGURES_DIR / "figure4_mapped_difference_pct_by_depth.png"
    ),
    "figure4_pdf": (
        FIGURES_DIR / "figure4_mapped_difference_pct_by_depth.pdf"
    ),
    "figure5_png": (
        FIGURES_DIR / "figure5_sample_mapped_reads_by_depth.png"
    ),
    "figure5_pdf": (
        FIGURES_DIR / "figure5_sample_mapped_reads_by_depth.pdf"
    ),
    "summary": FIGURES_DIR / "figure_summary_stats.txt",
}


def reads_in_millions(value: float, position: float) -> str:
    """Format read counts as whole millions on plot axes."""
    del position
    return f"{value / 1_000_000:.0f}"


def load_and_validate_metrics() -> pd.DataFrame:
    """Load and validate the canonical Experiment 2a metrics table."""
    if not METRICS_TABLE.is_file():
        raise FileNotFoundError(
            f"Metrics table not found: {METRICS_TABLE}\n"
            "Run scripts/generate_metrics_table.py first."
        )

    try:
        metrics = pd.read_csv(METRICS_TABLE, sep="\t")
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

    observed_samples = set(metrics["sample"].astype(str))

    if observed_samples != EXPECTED_SAMPLES:
        raise ValueError(
            "Unexpected sample set in metrics table.\n"
            f"Observed: {sorted(observed_samples)}\n"
            f"Expected: {sorted(EXPECTED_SAMPLES)}"
        )

    duplicate_mask = metrics.duplicated(
        subset=["sample", "depth_percent"],
        keep=False,
    )

    if duplicate_mask.any():
        duplicate_rows = metrics.loc[
            duplicate_mask,
            ["sample", "depth_percent"],
        ]

        raise ValueError(
            "Duplicate sample-depth combinations found:\n"
            f"{duplicate_rows.to_string(index=False)}"
        )

    for sample in sorted(EXPECTED_SAMPLES):
        observed_depths = set(
            metrics.loc[
                metrics["sample"] == sample,
                "depth_percent",
            ].astype(int)
        )

        if observed_depths != EXPECTED_DEPTHS:
            raise ValueError(
                f"Unexpected depth series for sample {sample}.\n"
                f"Observed: {sorted(observed_depths)}\n"
                f"Expected: {sorted(EXPECTED_DEPTHS)}"
            )

    numeric_columns = [
        "depth_fraction",
        "depth_percent",
        "flagstat_total",
        "flagstat_mapped",
        "idxstats_total",
        "idxstats_mapped",
        "mapped_diff",
        "mapped_diff_pct",
    ]

    for column in numeric_columns:
        if metrics[column].isna().any():
            raise ValueError(
                f"Column contains missing values: {column}"
            )

    return (
        metrics.sort_values(
            ["sample", "depth_percent"],
            kind="stable",
        )
        .reset_index(drop=True)
    )


def configure_plotting() -> None:
    """Apply the publication-style plotting configuration."""
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 13,
            "legend.fontsize": 11,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "figure.dpi": 300,
        }
    )


def save_figure(png_path: Path, pdf_path: Path) -> None:
    """Save the current figure in PNG and PDF formats."""
    plt.tight_layout()
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()

    for path in (png_path, pdf_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise OSError(
                f"Figure output was not created or is empty: {path}"
            )


def generate_figure2(
    summary: pd.DataFrame,
) -> None:
    """Generate Figure 2: mean mapped reads across depth."""
    plt.figure(figsize=(7, 5))

    plt.plot(
        summary["depth_percent"],
        summary["flagstat_mapped"],
        marker="o",
        linewidth=2.4,
        label="samtools flagstat",
    )

    plt.plot(
        summary["depth_percent"],
        summary["idxstats_mapped"],
        marker="s",
        linestyle="--",
        linewidth=2.0,
        label="samtools idxstats",
    )

    plt.xlabel("Downsampling depth (%)")
    plt.ylabel("Mean mapped read count (millions)")
    plt.title("Mapped-read counts: flagstat versus idxstats")
    plt.gca().yaxis.set_major_formatter(
        FuncFormatter(reads_in_millions)
    )
    plt.grid(True, alpha=0.3)
    plt.legend(frameon=True)

    save_figure(
        FIGURE_FILES["figure2_png"],
        FIGURE_FILES["figure2_pdf"],
    )


def generate_figure3(
    metrics: pd.DataFrame,
    pearson_r: float,
    mean_absolute_difference: float,
    comparison_count: int,
) -> None:
    """Generate Figure 3: flagstat–idxstats agreement."""
    plt.figure(figsize=(6.8, 6))

    for sample, subset in metrics.groupby(
        "sample",
        sort=True,
    ):
        plt.scatter(
            subset["flagstat_mapped"],
            subset["idxstats_mapped"],
            s=45,
            alpha=0.75,
            label=sample,
            edgecolors="black",
            linewidths=0.4,
        )

    minimum_value = min(
        metrics["flagstat_mapped"].min(),
        metrics["idxstats_mapped"].min(),
    )
    maximum_value = max(
        metrics["flagstat_mapped"].max(),
        metrics["idxstats_mapped"].max(),
    )

    plt.plot(
        [minimum_value, maximum_value],
        [minimum_value, maximum_value],
        linestyle="--",
        linewidth=1.5,
        label="Identity line",
    )

    plt.xlabel("Mapped read count (flagstat, millions)")
    plt.ylabel("Mapped read count (idxstats, millions)")
    plt.title(
        "Agreement between flagstat and idxstats mapped-read counts"
    )

    plt.gca().xaxis.set_major_formatter(
        FuncFormatter(reads_in_millions)
    )
    plt.gca().yaxis.set_major_formatter(
        FuncFormatter(reads_in_millions)
    )

    plt.grid(True, alpha=0.3)

    plt.text(
        0.05,
        0.95,
        (
            f"Pearson r = {pearson_r:.4f}\n"
            f"Mean absolute difference = "
            f"{mean_absolute_difference:.0f}\n"
            f"n = {comparison_count}"
        ),
        transform=plt.gca().transAxes,
        verticalalignment="top",
        bbox={
            "boxstyle": "round",
            "alpha": 0.15,
        },
    )

    handles, labels = plt.gca().get_legend_handles_labels()

    identity_index = labels.index("Identity line")
    legend_order = [identity_index] + [
        index
        for index, label in enumerate(labels)
        if label != "Identity line"
    ]

    plt.legend(
        [handles[index] for index in legend_order],
        [labels[index] for index in legend_order],
        frameon=True,
        loc="lower right",
    )

    save_figure(
        FIGURE_FILES["figure3_png"],
        FIGURE_FILES["figure3_pdf"],
    )


def generate_figure4(
    metrics: pd.DataFrame,
) -> None:
    """Generate Figure 4: mapped-read percentage difference."""
    markers = ["o", "s", "^"]

    plt.figure(figsize=(7, 5))

    for index, (sample, subset) in enumerate(
        metrics.groupby("sample", sort=True)
    ):
        plt.plot(
            subset["depth_percent"],
            subset["mapped_diff_pct"],
            marker=markers[index],
            linewidth=2,
            label=sample,
        )

    plt.axhline(
        0,
        linestyle="--",
        linewidth=1.5,
    )

    plt.xlabel("Downsampling depth (%)")
    plt.ylabel("Mapped-read difference (%)")
    plt.title(
        "Mapped-read percentage difference by downsampling depth"
    )
    plt.ylim(-0.002, 0.002)
    plt.grid(True, alpha=0.3)

    plt.text(
        0.03,
        0.08,
        "All 27 comparisons showed 0.0% difference",
        transform=plt.gca().transAxes,
        bbox={
            "boxstyle": "round",
            "alpha": 0.15,
        },
    )

    plt.legend(frameon=True)

    save_figure(
        FIGURE_FILES["figure4_png"],
        FIGURE_FILES["figure4_pdf"],
    )


def generate_figure5(
    metrics: pd.DataFrame,
) -> None:
    """Generate Figure 5: per-sample mapped reads across depth."""
    markers = ["o", "s", "^"]

    plt.figure(figsize=(8, 5))

    for index, (sample, subset) in enumerate(
        metrics.groupby("sample", sort=True)
    ):
        plt.plot(
            subset["depth_percent"],
            subset["flagstat_mapped"],
            marker=markers[index],
            linewidth=2,
            label=sample,
        )

    plt.xlabel("Downsampling depth (%)")
    plt.ylabel("Mapped read count (millions)")
    plt.title("Mapped read counts across downsampling depths")

    plt.gca().yaxis.set_major_formatter(
        FuncFormatter(reads_in_millions)
    )

    plt.grid(True, alpha=0.3)
    plt.legend(frameon=True)

    save_figure(
        FIGURE_FILES["figure5_png"],
        FIGURE_FILES["figure5_pdf"],
    )


def write_summary_statistics(
    pearson_r: float,
    mean_absolute_difference: float,
    maximum_absolute_difference: float,
    comparison_count: int,
) -> None:
    """Write the figure-level quantitative summary."""
    summary_path = FIGURE_FILES["summary"]

    with summary_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            f"Pearson correlation: {pearson_r:.6f}\n"
        )
        handle.write(
            "Mean absolute mapped-read difference: "
            f"{mean_absolute_difference:.6f}\n"
        )
        handle.write(
            "Maximum absolute mapped-read difference: "
            f"{maximum_absolute_difference}\n"
        )
        handle.write(
            f"Number of comparisons: {comparison_count}\n"
        )

    if not summary_path.is_file() or summary_path.stat().st_size == 0:
        raise OSError(
            f"Summary statistics output was not created: {summary_path}"
        )


def main() -> None:
    """Generate and verify all Experiment 2a figures."""
    metrics = load_and_validate_metrics()
    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    configure_plotting()

    pearson_r = metrics["flagstat_mapped"].corr(
        metrics["idxstats_mapped"]
    )
    maximum_absolute_difference = (
        metrics["mapped_diff"].abs().max()
    )
    mean_absolute_difference = (
        metrics["mapped_diff"].abs().mean()
    )
    comparison_count = len(metrics)

    summary = (
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

    generate_figure2(summary)

    generate_figure3(
        metrics=metrics,
        pearson_r=pearson_r,
        mean_absolute_difference=mean_absolute_difference,
        comparison_count=comparison_count,
    )

    generate_figure4(metrics)
    generate_figure5(metrics)

    write_summary_statistics(
        pearson_r=pearson_r,
        mean_absolute_difference=mean_absolute_difference,
        maximum_absolute_difference=maximum_absolute_difference,
        comparison_count=comparison_count,
    )

    print(f"Wrote publication figures to {FIGURES_DIR}")
    print(f"Pearson correlation: {pearson_r:.6f}")
    print(
        "Mean absolute mapped-read difference: "
        f"{mean_absolute_difference:.6f}"
    )
    print(
        "Maximum absolute mapped-read difference: "
        f"{maximum_absolute_difference}"
    )
    print(f"Number of comparisons: {comparison_count}")

    for path in FIGURE_FILES.values():
        print(f"Verified output: {path}")


if __name__ == "__main__":
    try:
        main()
    except (
        FileNotFoundError,
        OSError,
        ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc