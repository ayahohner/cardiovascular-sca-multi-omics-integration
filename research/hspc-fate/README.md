# HSPC Lineage-Priming Module

## Status

The v2 workflow is implemented. Its gated result remains `insufficient_evidence`, now for stronger reasons than v1: held-out endpoint markers do not validate every healthy trajectory, reasonable model refits change effect directions, the two estimators disagree in validation, and validation donor SCD5 retains one reference-supported HSC/MPP-like cell against the frozen minimum of 20.

## Question and claim boundary

Among transcriptomically HSC/MPP-like cells, is SCD associated with a reproducible donor-level difference in healthy-reference lineage priming toward myeloid, erythroid–megakaryocytic, or lymphoid programs?

The output is a cross-sectional association. It is not observed fate, temporal direction, lineage tracing, a causal SCD effect, or clinical prediction. The historical [`v1-spec.md`](v1-spec.md) records the original MVP contract; [`migration-v2.md`](migration-v2.md) explains why and how it changed.

## Implemented method

- Discovery remains `S-BSST1258`; validation remains `GSE133181`.
- Duplicate gene symbols are collapsed before scoring or PCA.
- Discovery uses a multi-gene HSPC admission rule and an ambient-profile audit; Scrublet call rates and score distributions are released.
- Construction and held-out validation marker panels are disjoint.
- Equal numbers of cells per healthy donor define each reference.
- Palantir supplies the primary cell-level trajectory estimate.
- CellRank 2 with a CytoTRACE kernel supplies a complementary potency-based sensitivity estimate.
- A leave-one-healthy-donor-out distance distribution defines out-of-reference rejection.
- SCD cells map into the frozen healthy reference and cannot refit it.
- Donors remain the biological replicates.
- Donor deletions and two bounded resolution specifications refit the trajectory instead of merely recomputing means.

## Artifacts

- [`study-report.md`](study-report.md) — complete v2 scientific report.
- [`migration-v2.md`](migration-v2.md) — MVP-to-step-7 implementation record.
- [`operations.md`](operations.md) — current executable sequence and stop conditions.
- [`analysis.json`](analysis.json) — frozen source, QC, marker-role, model, and sensitivity configuration.
- [`sources.json`](sources.json) — required sources plus open and controlled step-7 pointers.
- [`../../results/hspc-fate/report.md`](../../results/hspc-fate/report.md) — compact released result and neighboring machine-readable artifacts.

## Reproduce

```bash
make setup
make hspc-fetch
make hspc-preflight
make hspc-mvp
make hspc-verify
make test
```

Raw archives and rebuildable intermediates remain outside Git. The compact release is committed under `results/hspc-fate/`.
