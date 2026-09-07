#!/usr/bin/env python3

"""
Generate Figure 1: overview of the reproducible computational
benchmarking framework.

Outputs:
    results/figures/figure1_benchmarking_workflow.png
    results/figures/figure1_benchmarking_workflow.pdf

The workflow represents the two complementary study components:

1. Baseline concordance experiment
2. Stale-index positive-control experiment

The figure is generated programmatically so that the study-design
schematic is reproducible and version controlled alongside the
analytical workflow.
"""

from __future__ import annotations

import sys
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
except ImportError as exc:
    raise SystemExit(
        "ERROR: matplotlib is required but is not installed.\n"
        "Create the project environment from environment.yml."
    ) from exc


RESULTS_DIR = Path("results")
FIGURES_DIR = RESULTS_DIR / "figures"

FIGURE_PNG = (
    FIGURES_DIR
    / "figure1_benchmarking_workflow.png"
)

FIGURE_PDF = (
    FIGURES_DIR
    / "figure1_benchmarking_workflow.pdf"
)


def add_box(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    fontsize: float = 10.2,
    linewidth: float = 1.3,
) -> None:
    """Add a rounded workflow box centred at x, y."""
    box = FancyBboxPatch(
        (
            x - width / 2,
            y - height / 2,
        ),
        width,
        height,
        boxstyle="round,pad=0.015",
        linewidth=linewidth,
        facecolor="white",
        edgecolor="black",
    )

    ax.add_patch(box)

    ax.text(
        x,
        y,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        wrap=True,
    )


def add_arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
) -> None:
    """Add a directional arrow between workflow elements."""
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=1.25,
        color="black",
        shrinkA=2,
        shrinkB=2,
    )

    ax.add_patch(arrow)


def save_figure() -> None:
    """Save Figure 1 in PNG and PDF formats."""
    plt.savefig(
        FIGURE_PNG,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.18,
    )

    plt.savefig(
        FIGURE_PDF,
        bbox_inches="tight",
        pad_inches=0.18,
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


def generate_workflow() -> None:
    """Generate the revised two-component study workflow."""
    fig, ax = plt.subplots(
        figsize=(12, 9)
    )

    ax.set_xlim(
        0,
        12,
    )

    ax.set_ylim(
        0.55,
        10.4,
    )

    ax.axis("off")

    # --------------------------------------------------------------
    # Common source and validation stage
    # --------------------------------------------------------------

    add_box(
        ax,
        6.0,
        9.65,
        4.35,
        0.95,
        (
            "Three publicly available prostate tumour\n"
            "whole-exome sequencing (WES) BAMs\n"
            "SRR10353673 • SRR10353687 • SRR10353710"
        ),
    )

    add_box(
        ax,
        6.0,
        8.15,
        4.45,
        0.95,
        (
            "Source-data verification\n"
            "ENA MD5 checksum • samtools quickcheck\n"
            "BAM index verification"
        ),
    )

    add_arrow(
        ax,
        (6.0, 9.16),
        (6.0, 8.65),
    )

    # --------------------------------------------------------------
    # Branch headings
    # --------------------------------------------------------------

    ax.text(
        3.0,
        7.05,
        "Baseline concordance experiment",
        ha="center",
        va="center",
        fontsize=11.5,
        fontweight="bold",
    )

    ax.text(
        9.0,
        7.05,
        "Stale-index positive-control experiment",
        ha="center",
        va="center",
        fontsize=11.5,
        fontweight="bold",
    )

    # --------------------------------------------------------------
    # Branch geometry
    #
    # A central vertical stem descends from the source-verification
    # box to a horizontal fork positioned below the branch headings.
    # Short vertical arrows then enter the first box in each branch.
    #
    # This keeps all connector lines clear of the heading text.
    # --------------------------------------------------------------

    fork_y = 6.70

    # Central stem from the verification box to the fork.
    ax.plot(
        [
            6.0,
            6.0,
        ],
        [
            7.67,
            fork_y,
        ],
        color="black",
        linewidth=1.25,
        solid_capstyle="butt",
    )

    # Horizontal fork.
    ax.plot(
        [
            3.0,
            9.0,
        ],
        [
            fork_y,
            fork_y,
        ],
        color="black",
        linewidth=1.25,
        solid_capstyle="butt",
    )

    # Short downward branch arrows into the first boxes.
    add_arrow(
        ax,
        (3.0, fork_y),
        (3.0, 6.47),
    )

    add_arrow(
        ax,
        (9.0, fork_y),
        (9.0, 6.47),
    )

    # --------------------------------------------------------------
    # Baseline branch
    # --------------------------------------------------------------

    add_box(
        ax,
        3.0,
        5.95,
        4.25,
        1.0,
        (
            "Nine deterministic nested downsampling\n"
            "conditions per source BAM\n"
            "10% • 20% • … • 90%"
        ),
    )

    add_box(
        ax,
        3.0,
        4.45,
        4.25,
        0.95,
        (
            "27 paired baseline observations\n"
            "3 source BAMs × 9 nested depth conditions"
        ),
    )

    add_box(
        ax,
        3.0,
        3.0,
        4.25,
        1.0,
        (
            "Mapped-read extraction\n"
            "samtools flagstat → F\n"
            "samtools idxstats → I"
        ),
    )

    add_box(
        ax,
        3.0,
        1.45,
        4.25,
        1.05,
        (
            "Baseline discrepancy assessment\n"
            "Discrepancy measures: Δ = F − I, |Δ|, and %Δ\n"
            "Exact-concordance count"
        ),
        fontsize=9.9,
    )

    add_arrow(
        ax,
        (3.0, 5.44),
        (3.0, 4.94),
    )

    add_arrow(
        ax,
        (3.0, 3.97),
        (3.0, 3.51),
    )

    add_arrow(
        ax,
        (3.0, 2.49),
        (3.0, 1.99),
    )

    # --------------------------------------------------------------
    # Positive-control branch
    # --------------------------------------------------------------

    add_box(
        ax,
        9.0,
        5.95,
        4.25,
        1.0,
        (
            "Representative validated baseline BAM\n"
            "SRR10353673 at 10% depth"
        ),
    )

    add_box(
        ax,
        9.0,
        4.45,
        4.25,
        1.15,
        (
            "Prespecified stale-index perturbations\n"
            "Remove k = 1,000 or 10,000\n"
            "mapped alignment records\n"
            "Retain the original BAM index"
        ),
        fontsize=9.7,
    )

    add_box(
        ax,
        9.0,
        3.0,
        4.25,
        1.0,
        (
            "Mapped-read extraction\n"
            "samtools flagstat → F\n"
            "samtools idxstats → I"
        ),
    )

    add_box(
        ax,
        9.0,
        1.45,
        4.25,
        1.05,
        (
            "Positive-control recovery assessment\n"
            "Expected Δ compared with observed Δ\n"
            "Residual R = Δobserved − Δexpected"
        ),
        fontsize=9.8,
    )

    add_arrow(
        ax,
        (9.0, 5.44),
        (9.0, 5.04),
    )

    add_arrow(
        ax,
        (9.0, 3.87),
        (9.0, 3.51),
    )

    add_arrow(
        ax,
        (9.0, 2.49),
        (9.0, 1.99),
    )

    save_figure()


def main() -> None:
    """Generate and verify Figure 1."""
    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    generate_workflow()

    print(
        "Generated reproducible study-workflow figure "
        f"in {FIGURES_DIR}."
    )

    for path in (
        FIGURE_PNG,
        FIGURE_PDF,
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
    except OSError as exc:
        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc