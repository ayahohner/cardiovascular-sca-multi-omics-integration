# HSPC Fate Module

## Research question

Among human HSC/MPP-like cells, is SCA associated with a reproducible difference in inferred allocation toward myeloid, erythroid–megakaryocytic, or lymphoid progenitor fates?

The estimand is the SCA-minus-healthy difference between donor means of inferred terminal-fate probability. It is an association derived from cross-sectional transcriptomes. It is not observed fate, lineage tracing, a causal disease effect, or a clinical prediction.

## Valid-result contract

The discovery cohort is `S-BSST1258`; the independent human validation cohort is `GSE133181`. Healthy cells define the reference topology. SCA cells are mapped without refitting that topology. All comparisons use donors, not cells, as biological replicates.

The workflow may return exactly one conclusion:

1. `replicated_inferred_allocation_difference` — at least one predeclared lineage has concordant discovery and validation direction and passes stability controls;
2. `dataset_specific_inferred_allocation_difference` — a discovery contrast passes internal controls but does not reproduce;
3. `no_stable_inferred_allocation_difference` — estimable contrasts do not pass predeclared stability and replication controls; or
4. `insufficient_evidence` — a source, topology, annotation, or inference gate prevents an interpretable contrast.

Every conclusion is a valid result when its gates, provenance, uncertainty, and failures are reported. Statistical significance is not a release condition because donor counts are small.

## Boundary

In scope: open processed transcriptomic data, three lineage families, two independent fate estimators, donor-level contrasts, exact label permutations, leave-one-donor-out sensitivity, independent human replication, and optional orthogonal program checks.

Out of scope: the endothelial-transition MVP, joint atlases, RNA velocity without compatible layers, direct lineage tracing, controlled-access data, treatment response, gene therapy efficacy, cardiovascular tissue inference, causal claims, individual prediction, and clinical use.

## Module artifacts

- [`sources.json`](sources.json) records required, optional, and excluded source pointers.
- [`operations.md`](operations.md) is the implementation-order and stop-gate contract.
- [`v1-spec.md`](v1-spec.md) is the formal EARS and Gherkin specification.

Future implementation is confined to `data/raw/hspc-fate/`, `data/derived/hspc-fate/`, and `results/hspc-fate/`. No HSPC asset is currently downloaded or committed.
