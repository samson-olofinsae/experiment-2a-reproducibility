# Experiment 2a Reproducibility Package

## Overview

This repository contains the complete computational workflow used to reproduce **Experiment 2a** described in the manuscript:

> **A Reproducible Computational Benchmarking Framework for Evaluating Agreement Between BAM-Derived Mapped-Read Metrics**

Experiment 2a establishes a reproducible computational benchmarking workflow for evaluating agreement between mapped-read counts reported by **samtools flagstat** and **samtools idxstats** using three publicly available prostate tumour BAM datasets analysed across nine sequencing-depth fractions (10–90%).

This repository serves as the archived reference implementation of the computational workflow reported in Experiment 2a and complements the separate bam-check software repository, which provides the general-purpose implementation of the benchmarking methodology.

The repository enables reviewers and researchers to reproduce the complete computational workflow, including:

- downloading the three public BAM datasets from the European Nucleotide Archive (ENA);
- verifying dataset integrity using MD5 checksum validation;
- generating deterministic downsampled BAM datasets;
- extracting mapped-read metrics using **samtools flagstat** and **samtools idxstats**;
- generating the benchmarking tables;
- producing the publication figures;
- comparing reproduced outputs against the archived expected results.

---

# Repository Structure

```
experiment-2a-reproducibility/

config/
    experiment_2a_bam_manifest.tsv

scripts/
    download_bams.sh
    downsample.sh
    compute_flagstat.sh
    compute_idxstats.sh
    generate_metrics_table.py
    plot_metrics.py
    validate_pipeline.sh
    run_full_analysis.sh
    create_supplementary_datasets.py

expected_results/
    metrics_table.tsv
    discrepancy_summary.tsv
    figures/

supplementary/
    Supplementary_Dataset_S1.xlsx
    Supplementary_Dataset_S2.xlsx

environment.yml
LICENSE
README.md
```

---

# Software Requirements

The computational workflow was developed and tested using the Conda environment provided in:

```
environment.yml
```

Create the environment using:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate experiment-2a-reproducibility
```

The environment installs all software required to reproduce the workflow, including:

- Python 3.10
- samtools 1.22
- pandas
- matplotlib
- aria2
- openpyxl

---

# Input Datasets

Experiment 2a uses three publicly available prostate tumour BAM datasets obtained from the **European Nucleotide Archive (ENA)**.

| Sample | ENA Study |
|---------|-----------|
| SRR10353673 | PRJNA577801 |
| SRR10353687 | PRJNA577801 |
| SRR10353710 | PRJNA577801 |

The file

```
config/experiment_2a_bam_manifest.tsv
```

contains the download URLs and MD5 checksums used during the study.

---

# Running the Workflow

Download the BAM datasets:

```bash
bash scripts/download_bams.sh
```

Run the complete Experiment 2a workflow:

```bash
bash scripts/run_full_analysis.sh
```

The workflow performs:

1. dataset validation;
2. deterministic BAM downsampling;
3. mapped-read metric extraction;
4. benchmarking table generation;
5. publication figure generation.

---

# Expected Outputs

Successful execution produces:

```
results/

metrics_table.tsv

discrepancy_summary.tsv

figures/

    figure2_mean_mapped_reads_by_depth.pdf
    figure2_mean_mapped_reads_by_depth.png

    figure3_flagstat_vs_idxstats_agreement.pdf
    figure3_flagstat_vs_idxstats_agreement.png

    figure4_mapped_difference_pct_by_depth.pdf
    figure4_mapped_difference_pct_by_depth.png

    figure5_sample_mapped_reads_by_depth.pdf
    figure5_sample_mapped_reads_by_depth.png

    figure_summary_stats.txt
```

---

# Verifying Reproducibility

Reference outputs generated during the study are provided in:

```
expected_results/
```

Successful reproduction of Experiment 2a should produce benchmarking tables and publication figures consistent with these archived reference outputs.

Minor differences in PDF metadata (for example, creation timestamps) may occur depending on the local computing environment; however, the numerical benchmarking results should reproduce the published values.

---

# Supplementary Datasets

Excel versions of the principal benchmarking datasets are provided for convenience.

```
supplementary/

Supplementary_Dataset_S1.xlsx
Supplementary_Dataset_S2.xlsx
```

These correspond directly to:

- `metrics_table.tsv`
- `discrepancy_summary.tsv`

respectively.

---

# Citation

If you use this computational workflow in your research, please cite the accompanying manuscript:

> Olofinsae SA, Fatumo S, *et al.*
>
> **A Reproducible Computational Benchmarking Framework for Evaluating Agreement Between BAM-Derived Mapped-Read Metrics.**
>
> *(Publication details to be updated upon acceptance.)*

---

# License

This repository is distributed under the MIT License. See the accompanying `LICENSE` file for details.