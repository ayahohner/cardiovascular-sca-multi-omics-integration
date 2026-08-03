# Open-Data HSPC Lineage-Priming Study: v2 Report

- **Project:** Cardiovascular SCA Multi-Omics Integration
- **Module:** HSPC lineage priming
- **Release:** v2
- **Classification:** `insufficient_evidence`

## Executive summary

This study asks whether open human single-cell data support a reproducible association between sickle cell disease and HSC/MPP lineage priming toward myeloid, erythroid–megakaryocytic, or lymphoid programs.

Version 2 replaces the MVP's marker-derived centroid graph with a donor-balanced, cell-level healthy reference. Palantir is the primary trajectory estimator; CellRank 2 with a CytoTRACE kernel is a complementary potency-based sensitivity analysis. Construction markers are disjoint from markers used to evaluate the resulting root and endpoints. SCD cells map into the frozen healthy reference and may be rejected as out of distribution. Donor and bounded hyperparameter perturbations refit the model rather than only recomputing group means.

The result remains `insufficient_evidence`. This is not a failure to obtain a preferred answer. It is the defensible outcome of stronger tests:

1. the held-out erythroid–megakaryocytic endpoint fails in discovery;
2. held-out erythroid, lymphoid, and root checks fail in validation;
3. reasonable donor and resolution refits change effect directions or fail to produce a valid connected model;
4. Palantir and CellRank disagree on two of three validation directions; and
5. validation donor SCD5 retains one reference-supported HSC/MPP-like cell, below the frozen minimum of 20.

## Research question and estimand

The precise question is:

> Among transcriptomically HSC/MPP-like cells, is SCD associated with a reproducible donor-level difference in healthy-reference lineage priming?

For each donor, estimator, and lineage, the workflow averages lineage probabilities across that donor's reference-supported HSC/MPP-like cells. The condition effect is the mean SCD donor estimate minus the mean healthy donor estimate. Donors—not cells—are the biological replicates.

“Lineage priming” means position in a computationally directed healthy expression reference. It is not an observed descendant, longitudinal transition, causal disease effect, or calibrated clinical probability.

## Data

| Role | Accession | Comparison | Context |
|---|---|---|---|
| Discovery | [S-BSST1258](https://www.ebi.ac.uk/biostudies/studies/S-BSST1258) | 3 healthy and 3 SCD donors | Plerixafor-mobilized CD34-positive HSPCs in a gene-therapy collection context |
| Validation | [GSE133181](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133181) | 4 healthy and 5 SCD donors | Bone-marrow CD34-positive cells from a distinct pediatric study |

All seven required archives are open, checksum-verified, and donor resolved. Raw archives and rebuildable matrices remain outside Git. The release records URLs, byte counts, checksums, retrieval times, terms, assay, tissue, processing state, donor context, and limitations.

The collection contexts differ and are partly confounded with cohort. Cross-cohort transport is therefore a test, not an assumption.

## Method

### Preprocessing

Duplicate symbols are collapsed deterministically while retaining the first stable feature identifier. Discovery retains the source's 3,000-UMI floor, at least 500 genes, no more than 15% mitochondrial counts, and expression of at least one of `CD34`, `KIT`, or `ENG`. Validation retains its source-filtered-matrix rules. Scrublet is applied per donor, and observed call fractions and score distributions are released.

The discovery archives are unfiltered. Low-UMI barcodes define an ambient-expression profile used as a diagnostic, but exact EmptyDrops FDR is not reproduced because a maintained command-line-compatible implementation was not available in the current environment. This gap remains explicit.

The revised QC retained 64,751 discovery cells and 44,117 validation cells. These totals remain within the frozen source-comparison tolerance. Scrublet call fractions remain implausibly low for several discovery donors, so the doublet audit is a failed release gate.

### Healthy reference and biological anchors

Each healthy donor contributes 2,000 cells to the reference. PCA is fitted only on this balanced healthy subset using 1,000 variable genes and 25 components.

Five construction markers define the root anchor and each terminal family. Five disjoint markers per program evaluate the recovered root and endpoints. Held-out failure is preserved; markers are not swapped after viewing results.

### Trajectory estimators

Palantir uses a 30-neighbor diffusion graph, 10 diffusion components, and at most 600 waypoints. Three marker-anchored terminal seeds define the broad lineage families. Solver-level row-sum drift is permitted only below 2%, then probabilities are deterministically normalized and the observed deviation is released.

CellRank 2 constructs a separate CytoTRACE-directed transition kernel from the healthy reference and calculates probabilities toward donor-balanced terminal cell sets. Because both methods share input cells and biological endpoints, agreement is a sensitivity check rather than independent biological replication.

### Out-of-reference rejection

Each healthy donor is mapped against the other healthy donors. The 99th percentile of median 15-neighbor distances becomes the frozen OOD threshold. More than 98% of cells in every current donor are supported. Disease-query rejection is therefore implemented but is not the present cause of failure.

### Uncertainty and robustness

The workflow records exact donor-label permutations, descriptive donor bootstraps, ordinary donor deletions, and estimator agreement. It also adds model-level perturbations:

- deleting a healthy donor refits gene selection, PCA, trajectory, endpoints, OOD calibration, and query mapping;
- deleting an SCD donor preserves the healthy model and recomputes the contrast; and
- lower- and higher-resolution specifications vary seed, healthy reference size, HVG count, PCA dimension, and graph neighborhood.

Failed refits are explicit data. They are not removed from the stability decision.

## Results

### Root and reference support

Every discovery donor passes the 20-root-cell floor. Eight of nine validation donors pass. SCD5 retains 157 of 160 cells inside the healthy OOD boundary but only one satisfies the root definition. Its failure is therefore biological/annotation support, not broad OOD rejection.

### Baseline condition effects

Effects are SCD minus healthy probability points.

| Cohort | Lineage | Palantir | CellRank-CytoTRACE | Direction agreement |
|---|---|---:|---:|---|
| Discovery | Myeloid | −0.0369 | −0.0001 | Yes |
| Discovery | Erythroid–megakaryocytic | +0.0327 | +0.0013 | Yes |
| Discovery | Lymphoid | +0.0042 | −0.0012 | No |
| Validation | Myeloid | −0.2276 | −0.0025 | Yes |
| Validation | Erythroid–megakaryocytic | +0.2708 | −0.0034 | No |
| Validation | Lymphoid | −0.0432 | +0.0059 | No |

Baseline agreement is not sufficient. Across completed discovery refits, Palantir myeloid effects range from −0.0739 to +0.0824, and erythroid–megakaryocytic effects range from −0.0862 to +0.0784. CellRank effects are smaller but also cross zero. No lineage passes the full model-robustness contract.

Three validation robustness specifications fail technically: two healthy-donor deletion refits and the lower-resolution specification. These failures are themselves evidence that the validation geometry is not reliably transportable under the current donor set.

### Gates and classification

The principal failed gates are:

- `doublet_call_rate_audit`;
- `model_robustness_execution`;
- `discovery_topology`;
- `validation_topology`; and
- `validation_root_cell_reliability`.

The deterministic classification is therefore `insufficient_evidence`. The analysis does not support a replicated SCD-associated lineage-priming difference. It also does not prove that no biological difference exists.

## Ethics, fairness, and interpretation

SCD genotype is not treated as a proxy for ancestry, race, geography, environment, treatment access, or socioeconomic context. Missing donor metadata remain `not_reported`. Differences between mobilized peripheral HSPCs and pediatric bone marrow are treated as transport limitations rather than automatically removed batch effects.

No individual prediction, treatment recommendation, gene-therapy efficacy claim, causal inference, or cardiovascular-tissue inference is permitted.

## Step-7 evidence search

No third donor-resolved healthy-versus-SCD cohort was confirmed as openly downloadable without approval.

- `EGAD50000001522` contains 3 healthy and 4 SCA peripheral-blood HSPC donors and is the most aligned validation candidate, but EGA requires controlled access. It is excluded under the active open-only boundary.
- `GSE285943` contains 148 healthy circulating-HSPC donors and is an attractive external healthy reference. It has no SCD comparator and cannot validate a condition effect.

No admitted source provides compatible SCD-specific time points, lineage tracing, metabolic labeling, or spliced/unspliced layers across both cohorts. Actual temporal direction therefore remains a data-acquisition problem.

## Reproduction and artifacts

```bash
make setup
make hspc-fetch
make hspc-preflight
make hspc-mvp
make hspc-verify
make test
```

The compact release includes source and donor manifests, QC and gate ledgers, trajectory definitions, OOD mapping summaries, donor lineage-priming estimates, full robustness results, stability and replication tables, claim boundaries, and a checksum manifest. See [`migration-v2.md`](migration-v2.md) for the detailed MVP-to-v2 change record.

## Conclusion

The modernized workflow is a substantive improvement even though it does not produce a positive finding. It breaks marker circularity, removes arbitrary centroid topology, balances healthy donors, rejects unsupported disease queries, and includes model refitting in uncertainty. Those stronger controls show that the current open cohorts cannot support a robust SCD-associated HSPC lineage-priming result. The next useful increment is stronger data, not further tuning of this answer.
