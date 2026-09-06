# Experiment 2b: Positive-control discordance experiment

## Objective

To determine whether the benchmarking framework can quantitatively recover prespecified departures from samtools flagstat–samtools idxstats mapped-read concordance when BAM content and its corresponding index are deliberately rendered inconsistent.

## Software environment

- samtools: 1.22.1
- htslib: 1.23.1

## Experimental material

Representative baseline BAM:

SRR10353673_10pct.bam

Baseline mapped-read count:

M = 9,061,436

The positive-control derivatives are analytical controls and are not treated as independent biological datasets.

## Baseline condition

The unmodified BAM with its correct matching index must exhibit:

F = 9,061,436
I = 9,061,436
Delta = F - I = 0

This baseline condition was confirmed before positive-control generation.

## Perturbation mechanism

Mapped alignment records will be deterministically removed from copies of the baseline BAM while retaining the index generated from the unmodified baseline BAM.

Records will be removed in their original coordinate-order stream. Exactly the first k records satisfying the mapped-record criterion (SAM flag 0x4 absent) will be omitted. All remaining alignment records will be preserved in their original order.

The modified BAM will not be reindexed. Instead, the original baseline BAI will be associated deliberately with the modified BAM to create a controlled stale-index condition.

## Positive-control conditions

### PC1

k = 1,000 mapped alignment records removed.

Expected:

F = 9,060,436
I = 9,061,436
Delta_expected = -1,000

### PC2

k = 10,000 mapped alignment records removed.

Expected:

F = 9,051,436
I = 9,061,436
Delta_expected = -10,000

## Primary outcome

For each condition:

Delta_observed = F_observed - I_observed

Residual:

R = Delta_observed - Delta_expected

Exact quantitative recovery is defined as:

R = 0

## Success criterion

The baseline must remain concordant, and both stale-index positive controls must produce discordance in the prespecified direction and magnitude.

Expected states:

Baseline: Delta = 0
PC1: Delta = -1,000
PC2: Delta = -10,000

## Failure criterion

Any of the following will be retained and reported rather than replaced post hoc:

- persistence of concordance in a positive-control condition;
- incorrect discrepancy direction;
- incorrect discrepancy magnitude;
- unexpected samtools error;
- failure to remove exactly k mapped alignment records.

## Scope

This experiment tests one defined index-integrity failure mechanism at two prespecified perturbation magnitudes. It is not intended to characterize BAM pathology broadly, estimate diagnostic sensitivity or specificity, or replace the planned Synthetic BAM Benchmark Framework.
