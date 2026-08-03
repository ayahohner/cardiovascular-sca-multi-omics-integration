# HSPC Lineage-Priming v2 Operation Sequence

The donor is the unit of inference. A stopped or null result is valid when its provenance, estimates, uncertainty, and failed gates are released.

| Step | Operation | Required control | Failure response |
|---:|---|---|---|
| 1 | Freeze the estimand and marker roles | Construction and validation marker sets are disjoint; claims say healthy-reference lineage priming | Stop release until roles and language are unambiguous |
| 2 | Admit and preprocess donor matrices | Assets match checksums; duplicate symbols are collapsed; every exclusion and ambient/doublet diagnostic is recorded | Quarantine failed assets or classify `insufficient_evidence` |
| 3 | Fit the donor-balanced healthy trajectory | Palantir and CellRank-CytoTRACE produce finite normalized probabilities; held-out markers test root and endpoints | Retain auditable estimates and fail the topology gate |
| 4 | Map queries with abstention | SCD cells cannot alter model hashes; cells beyond the cross-donor healthy distance threshold are excluded | Report rejection and stop affected donor contrasts |
| 5 | Refit robustness checks | Every donor deletion and both bounded resolution specifications are recorded; positive directions must survive refitting | Mark the lineage unstable or fail robustness execution |
| 6 | Compare cohorts | Discovery and validation remain separate; donor means, exact permutations, deletion effects, and estimator agreement are reported | Preserve disagreement; never pool cohorts to rescue it |
| 7 | Acquire stronger evidence | Admit another open donor-resolved SCD cohort, or compatible temporal/kinetic layers, under the same provenance contract | Keep the current result and record the access or modality blocker |
| 8 | Release | Manifest checksums, gate ledger, claim boundary, and deterministic classification verify | Withhold an invalid release |

Classification remains one of `replicated_inferred_allocation_difference`, `dataset_specific_inferred_allocation_difference`, `no_stable_inferred_allocation_difference`, or `insufficient_evidence`. No interactive override or post-result threshold change is permitted.
