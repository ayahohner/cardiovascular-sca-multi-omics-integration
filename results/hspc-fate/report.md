# HSPC Fate v1 Result

## Classification

`insufficient_evidence`

This is a cross-sectional donor-level analysis of inferred lineage allocation. It does not observe cell fate and does not establish that SCA causes an allocation difference.

Discovery uses six open donor-separated unfiltered 10x matrices. Cell admission applies the source paper's fixed 3,000-UMI floor with a recorded barcode-rank knee diagnostic; it does not reproduce EmptyDrops FDR calculations. Scrublet replaces the paper's DoubletFinder implementation. The validation cohort uses its own frozen source-aware thresholds. These are explicit lean-MVP method deviations, not claims of an exact paper reproduction.

The primary estimator is absorption through a healthy-only directed state graph. The independent estimator projects cells to pseudotime-defined healthy terminal centroids without using the graph. Donor means are the biological replicates.

## Discovery

- myeloid: absorbing-Markov effect -0.025044; projection effect -0.004673; stable=false.
- erythroid_megakaryocytic: absorbing-Markov effect 0.022988; projection effect 0.008927; stable=false.
- lymphoid: absorbing-Markov effect 0.002055; projection effect -0.004254; stable=false.

No discovery lineage passed all frozen method-concordance, leave-one-donor-out, and exact-permutation null controls. `stability_results.tsv` preserves every donor-label assignment and every donor-specific leave-one-out effect for audit.

## Validation gate

Validation independently recovered the root and three terminal families. The donor reliability gate failed for SCD5 (4 root cells); estimates below the 20-root-cell floor are retained for audit but cannot support replication.

## Failed gates

- validation_root_cell_reliability

## Follow-up

Acquire an open independent human SCA HSPC cohort in which every comparator donor retains at least 20 transcriptomically HSC/MPP-like cells under the frozen marker semantics. Re-run the same donor-level estimand without pooling cohorts or changing endpoints. Optional bulk-human or mouse evidence may qualify programs but cannot rescue this human replication gate.
