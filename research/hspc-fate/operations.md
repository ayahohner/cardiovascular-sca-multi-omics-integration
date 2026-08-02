# HSPC Fate v1 Operation Sequence

## Result contract

The unit of inference is the donor. The primary output is the SCA-minus-healthy difference in mean inferred terminal-fate probability from HSC/MPP-like cells to each of three lineage families: myeloid, erythroid–megakaryocytic, and lymphoid. Discovery and validation are analyzed independently. A null or stopped result is valid when the relevant gates and uncertainty are reported.

## Discrete operations

| ID | Operation | Input | Output | Pass condition | Failure response |
|---|---|---|---|---|---|
| OP-001 | Admit sources | `sources.json` | frozen retrieval manifest | Both required accessions, terms, cohort labels, resolved assets, byte counts, and SHA-256 checksums are recorded | Return `insufficient_evidence` for an unavailable required source |
| OP-002 | Retrieve required matrices | retrieval manifest | immutable raw archives under `data/raw/hspc-fate/` | Downloaded bytes and checksums match the manifest | Quarantine the asset and stop dependent work |
| OP-003 | Build donor manifests | raw archives and repository metadata | one discovery table and one validation table | Every cell or matrix maps to one unique donor, condition, assay, tissue, and inclusion decision | Resolve the ambiguity or exclude the affected donor before expression analysis |
| OP-004 | Call cells and apply source-aware QC | discovery matrices and donor manifest | discovery expression object plus QC ledger | All six discovery donors remain, thresholds are frozen, and retained cells pass the fixed barcode-admission floor, library-size, mitochondrial, CD34, and doublet controls; the barcode-rank knee is recorded diagnostically | Revise only a failed technical rule, rerun all donors, and record the version; otherwise stop |
| OP-005 | Freeze the state vocabulary | healthy discovery cells and predeclared markers | root, intermediate, and terminal-state definitions | HSC/MPP root cells and all three terminal families are distinguishable in at least two healthy donors | Return `insufficient_evidence` for fate estimation |
| OP-006 | Fit the healthy discovery topology | QC-passing healthy discovery cells | frozen healthy manifold, pseudotime, and terminal-state model | The graph connects the root to all terminal families without condition data and passes marker-order and graph-connectivity controls | Change only predeclared graph sensitivity parameters or stop |
| OP-007 | Estimate discovery fate twice | frozen topology and all discovery donors | two cell-level terminal-probability matrices | A primary absorbing Markov estimator and an independent projection to pseudotime-defined terminal centroids each produce normalized probabilities for eligible HSC/MPP-like cells | Mark the failed lineage or method unresolved and stop its contrast |
| OP-008 | Aggregate discovery contrasts | eligible fate probabilities and donor manifest | donor-by-lineage estimates and SCA-minus-healthy contrasts | Every estimate is a donor mean and every lineage has 3 healthy and 3 SCD donor values | Stop the affected lineage; never substitute cells for donors |
| OP-009 | Test discovery stability | discovery donor estimates | every exact permutation assignment, every donor-specific leave-one-out effect, method concordance, and null controls | Direction agrees between estimators, at least 5 of 6 leave-one-donor-out effects retain direction, and permuted labels do not create a stronger systematic result | Classify the contrast as unstable |
| OP-010 | Rebuild the validation analysis | `GSE133181` matrices and frozen endpoint semantics | independent validation donor estimates | The healthy-only validation topology reaches the same three endpoint families and uses only the 4 healthy and 5 SCD donors | Return `insufficient_evidence` for replication |
| OP-011 | Compare cohorts | discovery and validation donor estimates | replication table | Each lineage records effect direction, magnitude, interval, method agreement, and cohort transport limits | Preserve disagreement; do not pool cohorts to erase it |
| OP-012 | Classify the result | all gates and replication table | one four-state conclusion | Exactly one of `replicated_inferred_allocation_difference`, `dataset_specific_inferred_allocation_difference`, `no_stable_inferred_allocation_difference`, or `insufficient_evidence` is selected by the decision table | Withhold release until classification is deterministic |
| OP-013 | Release an auditable result | classified result and all ledgers | versioned `results/hspc-fate/` package | Provenance resolves to both required accessions, donor-level tables and failures are present, and the narrative contains no observed-fate or causal claim | Withhold the affected artifact and release only after repair |

## Decision table

| Discovery | Independent validation | Conclusion |
|---|---|---|
| Passes stability; at least one lineage has concordant direction | Same lineage is estimable, directionally concordant, and stable | `replicated_inferred_allocation_difference` |
| Passes stability | Estimable but not concordant or stable | `dataset_specific_inferred_allocation_difference` |
| Estimable but no lineage passes stability | Estimable and no lineage passes replication | `no_stable_inferred_allocation_difference` |
| Any required source, annotation, topology, or inference gate prevents an interpretable contrast | Any | `insufficient_evidence` |

No p-value threshold changes this table. Exact permutations quantify how little the small donor cohorts can resolve; they do not manufacture precision.
