# HSPC Lineage-Priming v2 Result

## Classification

`insufficient_evidence`

This is a cross-sectional donor-level analysis of healthy-reference lineage priming. It does not observe cell fate and does not establish that SCA causes a transition change.

Discovery uses six open donor-separated unfiltered 10x matrices. Cell admission retains the source paper's fixed 3,000-UMI floor, adds an ambient-profile audit, uses a multi-gene HSPC identity rule instead of CD34 alone, collapses duplicate symbols, and records Scrublet score and call-rate diagnostics. Exact EmptyDrops FDR reproduction remains unavailable without R/DropletUtils.

The primary estimator is Palantir on a donor-balanced healthy cell reference. CellRank 2 with a CytoTRACE kernel is a complementary potency-based sensitivity estimator. Construction and held-out validation markers are disjoint. SCD cells cannot refit the reference and are rejected when their healthy-reference distance exceeds a leave-one-healthy-donor-out threshold. Donor means are the biological replicates.

## Discovery

- myeloid: Palantir effect -0.036907; CellRank-CytoTRACE effect -0.000104; stable=false.
- erythroid_megakaryocytic: Palantir effect 0.032671; CellRank-CytoTRACE effect 0.001346; stable=false.
- lymphoid: Palantir effect 0.004236; CellRank-CytoTRACE effect -0.001241; stable=false.

`stability_results.tsv` preserves every donor-label assignment and donor-specific deletion effect. `model_robustness.tsv` adds healthy-reference refits and bounded hyperparameter refits; a positive result must survive both.

## Validation gate

Validation independently recovered the root and three terminal families. The donor reliability gate failed for SCD5 (1 root cells); estimates below the 20-root-cell floor are retained for audit but cannot support replication.

## Failed gates

- doublet_call_rate_audit
- model_robustness_execution
- discovery_topology
- validation_topology
- validation_root_cell_reliability

## Follow-up

Verify a third donor-resolved open human SCD HSPC cohort or obtain compatible directional molecular layers. Re-run the frozen lineage-priming estimand without pooling cohorts or changing endpoints. Time-resolved, lineage-tracing, or compatible spliced/unspliced data are required before upgrading the claim toward temporal fate direction.
