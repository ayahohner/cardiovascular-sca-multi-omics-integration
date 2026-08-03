# From the HSPC MVP to the Step-7 Boundary

## Why v1 changed

The v1 release was reproducible and appropriately stopped, but its scientific model was too self-referential. The same marker programs contributed to root identity, maturity, pseudotime, terminal states, and endpoint checks. Thirty KMeans centroids imposed branch geometry, and the secondary projection reused the same PCA and marker-defined endpoints. Its neighbor-only sensitivity check passed even when reasonable seeds, HVG counts, PCA dimensions, or state counts changed effect signs.

v2 preserves the useful parts of v1—open-source provenance, healthy-only fitting, donor-level inference, independent cohorts, deterministic failure gates, and the no-causal-claim boundary—while replacing those weak layers.

## What was implemented

### Step 1: estimand and semantic separation

The active estimand is now **healthy-reference lineage priming**. Five construction markers per root or lineage select biological anchors. Five disjoint markers test the resulting root and endpoints. Cross-sectional output is never described as observed fate.

### Step 2: data and QC corrections

Stable feature identifiers are retained and duplicate symbols are collapsed before scoring. Discovery replaces the CD34-only expression filter with detection of at least one of `CD34`, `KIT`, or `ENG`. Low-UMI barcodes define an ambient-profile diagnostic. Scrublet expected rates, observed call fractions, score quantiles, warnings, and exclusions are released per donor.

An exact EmptyDrops reproduction was investigated but not added: the available Python package depends on a removed packaging interface, and R/DropletUtils is not installed. Pinning an obsolete Python stack would reduce reliability. The source 3,000-UMI floor remains explicit, and exact EmptyDrops is recorded as a remaining reproduction gap.

### Step 3: modern healthy trajectory

Each healthy donor contributes the same number of reference cells. PCA is fitted on that balanced reference. Palantir produces the primary cell-level pseudotime and lineage probabilities. CellRank 2 uses a CytoTRACE kernel and separately solves fate probabilities toward marker-anchored terminal sets. CellRank is a complementary sensitivity view, not falsely described as independent replication.

### Step 4: query abstention

The OOD threshold is the 99th percentile of median neighbor distances obtained by mapping each healthy donor only against the other healthy donors. SCD cells beyond that frozen distance are excluded. More than 98% of cells in every current donor are supported, so OOD rejection is not the present bottleneck.

### Step 5: refit uncertainty

The release now contains complete donor deletions. Deleting a healthy donor refits gene selection, PCA, trajectory, endpoints, OOD calibration, and query mapping; deleting an SCD donor reuses the frozen healthy reference and recomputes the contrast. Two additional specifications vary seed, reference size, HVG count, PCA dimension, and graph neighborhood. Failed refits are released as failures rather than omitted.

### Step 6: real-data result

The revised discovery effects are:

| Lineage | Palantir | CellRank-CytoTRACE | Baseline direction agreement |
|---|---:|---:|---|
| Myeloid | −0.0369 | −0.0001 | Yes |
| Erythroid–megakaryocytic | +0.0327 | +0.0013 | Yes |
| Lymphoid | +0.0042 | −0.0012 | No |

No lineage is stable after full refits. Across completed discovery perturbations, even the Palantir myeloid and erythroid–megakaryocytic effects span both signs. The discovery held-out erythroid endpoint fails. Validation fails held-out erythroid, lymphoid, and root checks; Palantir and CellRank disagree for two of three lineages; and SCD5 retains one root cell. The correct v2 classification is therefore `insufficient_evidence`.

### Step 7: acquisition boundary

The search found no confirmed third cohort that is both donor-resolved, healthy-versus-SCD, and open without approval. `EGAD50000001522` is highly aligned—3 healthy and 4 SCA peripheral-blood HSPC donors—but is controlled through EGA and is excluded. `GSE285943` is open and contains 148 healthy circulating-HSPC donors; it is a strong future external reference but cannot validate an SCD effect. Both pointers and their limitations are recorded in `sources.json`.

No admitted dataset currently provides a confirmed, compatible SCD-specific temporal, lineage-tracing, metabolic-labeling, or spliced/unspliced layer across both cohorts. Step 7 therefore remains data acquisition, not another model-tuning exercise.

## Net improvement

v2 does not manufacture a positive result. It makes failure more informative: marker circularity is broken, centroid topology is removed, disease queries can be rejected, model uncertainty is refitted, and the exact next data boundary is explicit. The released conclusion is more conservative and more defensible than v1.
