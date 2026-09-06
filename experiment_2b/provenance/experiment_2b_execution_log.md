# Experiment 2b execution log

## Experimental design

Experiment 2b evaluated whether the benchmarking framework could quantitatively recover prespecified departures from mapped-read concordance between samtools flagstat and samtools idxstats under a controlled stale-index condition.

The representative baseline BAM was:

SRR10353673_10pct.bam

Baseline mapped-read count:

M = 9,061,436

Two prespecified perturbation magnitudes were evaluated:

- PC1: removal of exactly 1,000 mapped alignment records
- PC2: removal of exactly 10,000 mapped alignment records

For each positive-control BAM, mapped records were removed deterministically from the original coordinate-order stream while the original baseline BAM index was deliberately retained.

## Baseline

Observed:

flagstat mapped = 9,061,436
idxstats mapped = 9,061,436
Delta_observed = 0

Expected:

Delta_expected = 0

Residual:

R = 0

## PC1

Mapped alignment records removed:

k = 1,000

The modified BAM passed samtools quickcheck.

Observed:

flagstat mapped = 9,060,436
idxstats mapped = 9,061,436
Delta_observed = -1,000

Expected:

Delta_expected = -1,000

Residual:

R = 0

PC1 therefore showed exact quantitative recovery of the prespecified discrepancy.

## PC2

Mapped alignment records removed:

k = 10,000

The modified BAM passed samtools quickcheck.

Observed:

flagstat mapped = 9,051,436
idxstats mapped = 9,061,436
Delta_observed = -10,000

Expected:

Delta_expected = -10,000

Residual:

R = 0

PC2 therefore showed exact quantitative recovery of the prespecified discrepancy.

## Infrastructure interruption during PC2 generation

The first attempt to generate PC2 was interrupted by a WSL filesystem failure. The filesystem became read-only and basic system commands subsequently returned input/output errors.

The partially written PC2 BAM failed samtools quickcheck because the BAM was truncated and lacked the expected EOF block. Its mapped-read counts were therefore treated as invalid and were not interpreted as experimental results.

The WSL instance was shut down and the host computer restarted. After restart, filesystem writability was confirmed and the invalid PC2 BAM was deleted.

PC2 was then rerun without any change to the preregistered design, perturbation magnitude, expected values, or analysis criteria.

## Final outcome

All three prespecified states were recovered exactly:

- Baseline: Delta = 0
- PC1: Delta = -1,000
- PC2: Delta = -10,000

For both positive controls:

Residual R = 0

Experiment 2b therefore demonstrated exact quantitative recovery of the prespecified stale-index discrepancies under the tested conditions.

## Scope

This experiment evaluated one defined index-integrity failure mechanism at two prespecified perturbation magnitudes. It does not constitute broad BAM-pathology characterization or an estimate of diagnostic sensitivity or specificity.
