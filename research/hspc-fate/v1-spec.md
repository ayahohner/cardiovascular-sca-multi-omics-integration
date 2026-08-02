# Feature 1: HSPC Fate Result Contract

The v1 result is a donor-level comparison of inferred lineage allocation in two independent public human cohorts. Cross-sectional transcriptomes support association and hypothesis generation, not observed fate or causation.

## Functional Requirements

- **[REQ-001]** The HSPC fate workflow shall estimate the SCA-minus-healthy difference between donor means of terminal-fate probability from HSC/MPP-like cells to myeloid, erythroid–megakaryocytic, and lymphoid lineages.
- **[REQ-002]** The HSPC fate workflow shall use the donor as the biological replicate for every condition contrast and uncertainty calculation.
- **[REQ-003]** The HSPC fate workflow shall classify each completed analysis as replicated_inferred_allocation_difference, dataset_specific_inferred_allocation_difference, no_stable_inferred_allocation_difference, or insufficient_evidence.
- **[REQ-004]** The HSPC fate workflow shall describe every estimate as a cross-sectional inferred-allocation association and exclude observed-fate, lineage-tracing, causal, individual, treatment-response, and clinical claims.

## Acceptance Criteria

- **[AC-001]** (References REQ-001)
  - **Given** eligible donor-level fate probabilities
  - **When** the primary contrast table is generated
  - **Then** it contains one SCA-minus-healthy donor-mean difference for each of the three named lineages

- **[AC-002]** (References REQ-002)
  - **Given** any condition contrast or uncertainty result
  - **When** its calculation record is inspected
  - **Then** its independent observations are donors and no cell is counted as an independent replicate

- **[AC-003]** (References REQ-003)
  - **Given** a completed gate ledger
  - **When** result classification runs
  - **Then** exactly one permitted classification is recorded

- **[AC-004]** (References REQ-004)
  - **Given** the released tables, figures, and narrative
  - **When** the claim-boundary check runs
  - **Then** inferred allocation is labeled as association and none of the excluded claim types is present

# Feature 2: Module and Source Boundary

Required sources are open human HSPC cohorts with donor identifiers. Optional sources may qualify programs but cannot rescue a failed human fate result.

## Functional Requirements

- **[REQ-005]** The HSPC fate workflow shall operate independently of the endothelial module's data, code, state vocabulary, configuration, and results.
- **[REQ-006]** The HSPC fate workflow shall use S-BSST1258 as the discovery cohort and GSE133181 as the independent validation cohort.
- **[REQ-007]** WHEN a source is admitted, the HSPC fate workflow shall record its accession, resolved asset URL, byte count, SHA-256 checksum, retrieval time, terms, donor labels, assay, tissue, preprocessing state, and limitations.
- **[REQ-008]** IF a required source lacks verified open access, integrity, donor labels, condition labels, or reuse terms, THEN the HSPC fate workflow shall stop its dependent operations and classify the result as insufficient_evidence.
- **[REQ-009]** WHERE optional orthogonal validation is included, the HSPC fate workflow shall restrict GSE277682 to human bulk-program direction and GSE277681 to within-age mouse-model support.
- **[REQ-010]** IF a source requires controlled access, lacks an explicit usable license, or lacks a comparator required by the estimand, THEN the HSPC fate workflow shall exclude that source from quantitative v1 analysis.

## Acceptance Criteria

- **[AC-005]** (References REQ-005)
  - **Given** a clean checkout containing both modules
  - **When** the HSPC fate dependency and path manifest is inspected
  - **Then** every read and write path is under the HSPC namespace or its own source ledger and no endothelial artifact is a dependency

- **[AC-006]** (References REQ-006)
  - **Given** the frozen v1 source manifest
  - **When** required cohort roles are queried
  - **Then** S-BSST1258 is the sole discovery cohort and GSE133181 is the sole validation cohort

- **[AC-007]** (References REQ-007)
  - **Given** an admitted source
  - **When** its frozen manifest record is validated
  - **Then** all eleven specified provenance fields are populated and the downloaded bytes match the recorded checksum

- **[AC-008]** (References REQ-008)
  - **Given** a required source with a failed access, integrity, label, or terms check
  - **When** the source gate runs
  - **Then** no dependent expression operation starts and the result is insufficient_evidence

- **[AC-009]** (References REQ-009)
  - **Given** optional GSE277682 or GSE277681 evidence
  - **When** its result record is inspected
  - **Then** its evidence role matches the permitted source-specific role and it does not alter human fate classification

- **[AC-010]** (References REQ-010)
  - **Given** a proposed source with any listed exclusion condition
  - **When** source admission runs
  - **Then** the source is marked excluded and contributes no quantitative input

# Feature 3: Donors, Cells, and Healthy Topology

Healthy cells define the geometry and endpoint semantics. Condition data cannot redefine the reference after contrasts are visible.

## Functional Requirements

- **[REQ-011]** WHEN donor manifests are built, the HSPC fate workflow shall assign every included matrix and cell to one unique donor, condition, assay, tissue, cohort, and inclusion decision.
- **[REQ-012]** WHEN discovery matrices are processed, the HSPC fate workflow shall apply frozen source-aware empty-droplet, library-size, mitochondrial, CD34-identity, and doublet controls and record donor-level retention.
- **[REQ-013]** WHILE the state vocabulary is being defined, the HSPC fate workflow shall use only healthy discovery cells and freeze HSC/MPP root, intermediate, and three terminal-lineage definitions before estimating condition contrasts.
- **[REQ-014]** IF HSC/MPP root cells or any terminal lineage cannot be distinguished in at least two healthy discovery donors, THEN the HSPC fate workflow shall stop fate estimation and classify the result as insufficient_evidence.
- **[REQ-015]** WHEN the discovery topology is fitted, the HSPC fate workflow shall use only QC-passing healthy discovery cells and verify root-to-terminal connectivity, marker order, donor support, and graph stability.

## Acceptance Criteria

- **[AC-011]** (References REQ-011)
  - **Given** an included discovery or validation observation
  - **When** its donor-manifest row is queried
  - **Then** exactly one value is present for each specified identity and inclusion field

- **[AC-012]** (References REQ-012)
  - **Given** the six discovery donor matrices
  - **When** preprocessing completes
  - **Then** the QC ledger records the frozen threshold, pass count, and exclusion count for every control and donor

- **[AC-013]** (References REQ-013)
  - **Given** an unfrozen state vocabulary
  - **When** root and lineage definitions are created
  - **Then** only healthy discovery cells contribute and the definitions are versioned before any SCA contrast is computed

- **[AC-014]** (References REQ-014)
  - **Given** healthy discovery annotations
  - **When** the state-coverage gate runs
  - **Then** each required root or terminal class has qualifying cells in at least two donors or the result is insufficient_evidence

- **[AC-015]** (References REQ-015)
  - **Given** a fitted healthy discovery topology
  - **When** topology validation runs
  - **Then** every terminal lineage is connected to the root and all four specified controls have recorded passing results

# Feature 4: Fate Estimation and Discovery Stability

Two estimators and donor perturbations test whether the inferred direction is an algorithmic or single-donor artifact.

## Functional Requirements

- **[REQ-016]** WHEN the healthy discovery topology passes validation, the HSPC fate workflow shall estimate normalized terminal-fate probabilities with a primary absorbing Markov estimator and an independently configured pseudotime-fate estimator.
- **[REQ-017]** WHEN SCA cells are mapped, the HSPC fate workflow shall preserve the frozen healthy transformation, topology, root definition, terminal definitions, and model parameters.
- **[REQ-018]** WHEN cell-level fate probabilities are aggregated, the HSPC fate workflow shall calculate one mean probability per donor and lineage before computing a condition contrast.
- **[REQ-019]** The HSPC fate workflow shall compute exact donor-label permutations, six leave-one-donor-out discovery fits, estimator-direction concordance, and label-permutation null controls for every estimable lineage.
- **[REQ-020]** IF estimator directions disagree or fewer than five of six leave-one-donor-out discovery fits retain the full-data direction, THEN the HSPC fate workflow shall mark the discovery lineage unstable.

## Acceptance Criteria

- **[AC-016]** (References REQ-016)
  - **Given** a passing healthy topology
  - **When** both fate estimators complete
  - **Then** each eligible HSC/MPP-like cell has three nonnegative probabilities that sum to one within numerical tolerance for each estimator

- **[AC-017]** (References REQ-017)
  - **Given** the frozen healthy model and SCA cells
  - **When** SCA mapping completes
  - **Then** hashes of all five frozen model components equal their pre-mapping hashes

- **[AC-018]** (References REQ-018)
  - **Given** eligible cell-level probabilities and donor labels
  - **When** the contrast input table is generated
  - **Then** it contains exactly one value per donor, estimator, and lineage and contains no cell-level replicate rows

- **[AC-019]** (References REQ-019)
  - **Given** an estimable discovery lineage
  - **When** stability analysis completes
  - **Then** all exact label assignments, six leave-one-donor-out fits, both estimator directions, and the null-control result are recorded

- **[AC-020]** (References REQ-020)
  - **Given** a discovery lineage with estimator disagreement or fewer than five direction-retaining leave-one-donor-out fits
  - **When** discovery status is assigned
  - **Then** the lineage is marked unstable and cannot support replicated_inferred_allocation_difference or dataset_specific_inferred_allocation_difference

# Feature 5: Independent Validation and Classification

Validation repeats endpoint semantics without pooling cohorts or tuning to discovery results.

## Functional Requirements

- **[REQ-021]** WHEN discovery stability analysis completes, the HSPC fate workflow shall rebuild an independent healthy-only topology for the four healthy and five SCD GSE133181 donors using the frozen lineage semantics.
- **[REQ-022]** IF the validation cohort cannot recover the HSC/MPP root and all three terminal lineages under the frozen semantics, THEN the HSPC fate workflow shall classify the result as insufficient_evidence.
- **[REQ-023]** WHEN discovery and validation estimates are compared, the HSPC fate workflow shall keep cohorts separate and report lineage-specific direction, magnitude, uncertainty, stability, estimator agreement, and transport limitations.
- **[REQ-024]** WHEN all required gates are evaluated, the HSPC fate workflow shall apply the operations decision table without manual override to select the result classification.

## Acceptance Criteria

- **[AC-021]** (References REQ-021)
  - **Given** completed discovery stability results and admitted GSE133181 matrices
  - **When** validation inference starts
  - **Then** the validation model uses only the four healthy donors for topology fitting, retains the five SCD donors for mapping, and uses the frozen lineage semantics

- **[AC-022]** (References REQ-022)
  - **Given** a validation state-coverage failure
  - **When** the validation gate runs
  - **Then** no replication claim is generated and the result is insufficient_evidence

- **[AC-023]** (References REQ-023)
  - **Given** estimable discovery and validation contrasts
  - **When** the replication table is generated
  - **Then** every lineage has separate cohort rows containing all six specified comparison fields

- **[AC-024]** (References REQ-024)
  - **Given** a complete gate ledger
  - **When** the classifier is run repeatedly
  - **Then** each run returns the same single decision-table classification without an interactive choice

# Feature 6: Uncertainty, Ethics, and Release

The release preserves weak, null, contradictory, and stopped evidence. Small cohorts and collection context remain visible at the point of interpretation.

## Functional Requirements

- **[REQ-025]** The HSPC fate workflow shall report donor counts, effect magnitudes, exact permutation results, leave-one-donor-out directions, estimator agreement, missing metadata, cohort context, and transport limitations for every lineage.
- **[REQ-026]** The HSPC fate workflow shall distinguish SCD genotype from ancestry, geography, environment, treatment context, access to care, and self-identified race and report unavailable attributes as missing.
- **[REQ-027]** WHEN a release candidate is prepared, the HSPC fate workflow shall publish a source manifest, donor manifest, QC ledger, state definitions, topology controls, donor-level fate table, stability results, replication table, gate ledger, result classification, and claim-boundary statement.
- **[REQ-028]** IF a release artifact lacks provenance, exposes donor-level sensitive metadata beyond the public source, omits a failed gate, or violates the claim boundary, THEN the HSPC fate workflow shall withhold the release until the defect is corrected and revalidated.

## Acceptance Criteria

- **[AC-025]** (References REQ-025)
  - **Given** any released lineage result
  - **When** its uncertainty record is inspected
  - **Then** all eight specified evidence fields are populated or explicitly marked unavailable

- **[AC-026]** (References REQ-026)
  - **Given** the cohort metadata and released narrative
  - **When** the ethics checklist runs
  - **Then** SCD genotype is not used as a proxy for any listed attribute and each unavailable attribute is marked missing

- **[AC-027]** (References REQ-027)
  - **Given** a release candidate
  - **When** its manifest is validated
  - **Then** all eleven specified artifacts are present, versioned, and linked

- **[AC-028]** (References REQ-028)
  - **Given** a release candidate with any listed defect
  - **When** the release gate runs
  - **Then** no affected artifact is published until the corrected candidate passes provenance, privacy, failure-disclosure, and claim-boundary checks
