# Reproducibility Package for BAM-Derived Mapped-Read Concordance Benchmarking

## Overview

This repository contains the computational workflow and reference outputs accompanying the manuscript:

> **A Reproducible Computational Benchmarking Framework for Evaluating Concordance Between BAM-Derived Mapped-Read Counts**

The study evaluates concordance between corresponding mapped-read quantities derived using **samtools flagstat** and **samtools idxstats** under defined analytical conditions.

The computational study comprises two complementary components:

1. **Baseline concordance experiment (Experiment 2a):** three independently sequenced publicly available prostate tumour whole-exome sequencing (WES) BAMs are evaluated across nine deterministic nested downsampling fractions (10–90%) per source BAM, yielding 27 paired mapped-read observations.

2. **Stale-index positive-control experiment (Experiment 2b):** a representative baseline BAM is subjected to prespecified perturbations in which exactly 1,000 or 10,000 mapped alignment records are removed while retaining the original BAM index. These controls evaluate whether known departures from mapped-read concordance are detected and quantitatively recovered.

The repository preserves the scripts, dataset metadata, reference outputs, positive-control truth definitions, and provenance records required to trace the reported computational results.

---

# Repository Structure

```text
experiment-2a-reproducibility/

config/
    experiment_2a_bam_manifest.tsv

scripts/
    download_bams.sh
    downsample.sh
    compute_flagstat.sh
    compute_idxstats.sh
    generate_metrics_table.py
    plot_workflow.py
    plot_metrics.py
    validate_pipeline.sh
    run_full_analysis.sh
    create_supplementary_datasets.py

expected_results/
    metrics_table.tsv
    discrepancy_summary.tsv
    figures/
        figure1_benchmarking_workflow.png
        figure1_benchmarking_workflow.pdf
        figure2_mean_mapped_reads_by_depth.png
        figure2_mean_mapped_reads_by_depth.pdf
        figure3_signed_discrepancy_by_depth.png
        figure3_signed_discrepancy_by_depth.pdf
        figure4_stale_index_expected_vs_observed_discrepancy.png
        figure4_stale_index_expected_vs_observed_discrepancy.pdf

experiment_2b/
    experiment_2b_prespecification.md
    positive_control_truth.tsv
    provenance/
        experiment_2b_execution_log.md
        experiment_2b_software_versions.txt
    results/
        experiment_2b_results.tsv
    scripts/
        remove_first_k_mapped.sh
        run_experiment_2b.sh
        plot_experiment_2b.py

supplementary/
    Supplementary_Dataset_S1.xlsx
    Supplementary_Dataset_S2.xlsx

environment.yml
LICENSE
README.md
```

---

# Software Requirements

The baseline computational workflow is defined by the Conda environment provided in:

```text
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

The released reproducibility environment includes:

- Python 3.10
- samtools 1.22
- pandas
- Matplotlib
- aria2
- openpyxl

Exact software versions used for the released analysis are documented in the reproducibility records. The stale-index positive-control workflow additionally uses GNU Awk.

---

# Input Datasets

The baseline concordance experiment uses three publicly available prostate tumour whole-exome sequencing BAM datasets obtained from the **European Nucleotide Archive (ENA)**.

| Sample | ENA Study | Library strategy |
|---|---|---|
| SRR10353673 | PRJNA577801 | WXS |
| SRR10353687 | PRJNA577801 | WXS |
| SRR10353710 | PRJNA577801 | WXS |

The manifest:

```text
config/experiment_2a_bam_manifest.tsv
```

records the source URLs and ENA reference MD5 checksums used for dataset acquisition and integrity verification.

---

# Baseline Concordance Experiment (Experiment 2a)

Each validated source BAM is downsampled directly to fractions representing 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, and 90% of the source alignment data using **samtools view -s** with a fixed sampling seed.

The fractions are generated directly from their respective source BAMs rather than sequentially from preceding fractions. Because the same deterministic sampling configuration is used across increasing fractions, the resulting alignment sets are nested. The nine fractions from each source BAM are therefore treated as derived analytical conditions rather than independent experimental replicates.

Mapped-read counts are then extracted using **samtools flagstat** and **samtools idxstats**. For the benchmark:

- **F** is the sum of QC-passed and QC-failed mapped records reported by `samtools flagstat`;
- **I** is the sum of mapped alignment counts reported across reference sequences by `samtools idxstats`.

The two utilities provide complementary alignment summaries; the benchmark compares only their corresponding mapped-read quantities and does not assume general interchangeability of the utilities.

---

# Concordance and Discrepancy Measures

For each paired observation, signed discrepancy is defined as:

```text
Δ = F - I
```

Absolute discrepancy is:

```text
|Δ| = |F - I|
```

and percentage discrepancy is:

```text
%Δ = |F - I| / F × 100
```

for observations where `F > 0`.

Baseline concordance is summarised using exact equality, absolute discrepancy, percentage discrepancy, maximum absolute discrepancy, and maximum percentage discrepancy.

The nine downsampling fractions within each source BAM are deterministic nested conditions and are not treated as independent replicates for population-level statistical inference. Agreement is assessed directly using discrepancy measures rather than Pearson correlation, because correlation does not establish equality or agreement between paired measurements.

---

# Running the Baseline Workflow

Download and validate the source BAM datasets:

```bash
bash scripts/download_bams.sh
```

Run the complete baseline analytical workflow:

```bash
bash scripts/run_full_analysis.sh
```

The baseline workflow performs:

1. source-data validation;
2. deterministic nested BAM downsampling;
3. mapped-read metric extraction;
4. benchmarking-table generation;
5. discrepancy analysis;
6. generation of the data-derived baseline publication figures.

Reference outputs are maintained in:

```text
expected_results/
```

The baseline workflow generates Figures 2 and 3 through:

```text
scripts/plot_metrics.py
```

Figure 1 is a reproducible study-design schematic generated separately through:

```bash
python scripts/plot_workflow.py
```

---

# Stale-Index Positive-Control Experiment (Experiment 2b)

The positive-control experiment evaluates whether the benchmarking framework can detect and quantitatively recover prespecified departures from mapped-read concordance under a defined stale-index condition.

The 10% downsampled BAM derived from **SRR10353673** serves as the representative baseline. Two derivative BAMs are constructed by removing exactly:

- 1,000 mapped alignment records; or
- 10,000 mapped alignment records

while retaining a copy of the original baseline BAM index.

The expected signed discrepancies are therefore prespecified as:

```text
PC1: Δexpected = -1,000
PC2: Δexpected = -10,000
```

Quantitative recovery is evaluated using:

```text
R = Δobserved - Δexpected
```

with exact recovery defined as `R = 0`.

The experimental design, expected values, success criteria, scripts, execution provenance, software versions, and observed results are preserved under:

```text
experiment_2b/
```

The positive controls evaluate the specified stale-index mechanism only and are not intended to establish general sensitivity or specificity across BAM pathologies.

To reproduce the stale-index positive-control experiment, use:

```bash
bash experiment_2b/scripts/run_experiment_2b.sh
```

Figure 4 is generated from the positive-control results using:

```bash
python experiment_2b/scripts/plot_experiment_2b.py
```

---

# Publication Figures

The four canonical publication figures are preserved in:

```text
expected_results/figures/
```

They are:

1. **Figure 1 — Benchmarking workflow:** generated by `scripts/plot_workflow.py`.
2. **Figure 2 — Mean mapped-read counts across nested depth conditions:** generated by `scripts/plot_metrics.py`.
3. **Figure 3 — Observation-level signed-discrepancy matrix:** generated by `scripts/plot_metrics.py`.
4. **Figure 4 — Expected versus observed stale-index discrepancy recovery:** generated by `experiment_2b/scripts/plot_experiment_2b.py`.

Figures 2 and 3 are data-derived outputs of the baseline analytical workflow. Figure 1 is a reproducible study-design schematic, and Figure 4 is generated from the separately reproduced stale-index positive-control experiment.

---

# Verifying Reproducibility

Successful reproduction should recover the archived numerical benchmarking outputs from the defined source data and computational environment.

For the baseline experiment, reproduction should recover the mapped-read counts and discrepancy measures stored in the reference analytical outputs.

For the positive-control experiment, the observed discrepancies should recover the prespecified truth values documented in:

```text
experiment_2b/positive_control_truth.tsv
```

and:

```text
experiment_2b/results/experiment_2b_results.tsv
```

Minor differences in non-analytical file metadata, such as PDF creation timestamps, may occur across computing environments; numerical results should remain consistent with the archived reference outputs.

---

# Supplementary Datasets

Machine-readable benchmarking outputs are provided in the repository, with Excel versions of the principal baseline datasets supplied for convenience:

```text
supplementary/
    Supplementary_Dataset_S1.xlsx
    Supplementary_Dataset_S2.xlsx
```

**Supplementary Dataset S1** contains the 27 paired baseline mapped-read observations and associated discrepancy measures.

**Supplementary Dataset S2** contains depth-level and overall summary statistics derived from the baseline benchmarking table.

These datasets are generated programmatically from the corresponding analytical outputs.

---

# Reproducibility and Provenance

The repository is maintained under Git version control. Scripts, analytical outputs, dataset metadata, positive-control truth definitions, and provenance records are retained so that reported results can be traced to their computational source.

The positive-control prespecification was created before execution of the corresponding positive-control experiments and records the perturbation mechanism, expected discrepancies, analysis criteria, and success criteria.

A versioned archival release of the repository is maintained through Zenodo. The release corresponding to the revised manuscript will be identified here after completion of the final reproducibility audit.

---

# Citation

If you use this computational workflow, please cite the accompanying manuscript:

> Olofinsae SA, Abubakar SD, Ojo I, Ojewumi OO, Corpas M, Fatumo S.
>
> **A Reproducible Computational Benchmarking Framework for Evaluating Concordance Between BAM-Derived Mapped-Read Counts**
>
> *(Publication details to be updated.)*

---

# License

This repository is distributed under the MIT License. See `LICENSE` for details.
