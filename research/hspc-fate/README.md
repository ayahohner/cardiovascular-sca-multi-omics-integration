# HSPC Fate Module

## Status

v1 is implemented and released. Its valid gated conclusion is `insufficient_evidence`: none of the three discovery contrasts passed all frozen stability controls, and validation donor SCD5 retained 4 HSC/MPP-like cells against the predeclared minimum of 20. Validation estimates remain available for audit but cannot support replication.

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
- [`analysis.json`](analysis.json) freezes required asset sizes and checksums, QC rules, lineage markers, model parameters, and claim limits.
- [`operations.md`](operations.md) is the implementation-order and stop-gate contract.
- [`v1-spec.md`](v1-spec.md) is the formal EARS and Gherkin specification.
- [`../../results/hspc-fate/report.md`](../../results/hspc-fate/report.md) reports the v1 result; its neighboring manifest and ledgers make the result auditable.

## Implemented method

Each cohort is analyzed independently. Healthy donors define a root, pseudotime, and three terminal-lineage families. The primary estimator uses absorption through the healthy-only state graph; the second projects to pseudotime-defined terminal centroids without using that graph. SCA cells are mapped through frozen model components, and only donor means enter contrasts, exact permutations, leave-one-donor-out checks, and bootstrap intervals.

This lean implementation does not claim an exact reproduction of either source paper. Discovery cell admission uses the published 3,000-UMI floor with a recorded barcode-rank knee diagnostic rather than rerunning EmptyDrops FDR, and Scrublet replaces DoubletFinder. These deviations, every permutation assignment, every leave-one-donor-out effect, QC retention, failed gates, and model-integrity hashes are released explicitly.

## Reproduce

```bash
make setup
make hspc-fetch
make hspc-preflight
make hspc-mvp
make hspc-verify
make test
```

Raw archives and rebuildable donor-level intermediates stay in `data/raw/hspc-fate/` and `data/derived/hspc-fate/` and remain outside Git. The compact release is committed under `results/hspc-fate/`.
