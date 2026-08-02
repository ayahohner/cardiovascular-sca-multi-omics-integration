# Feature 1: Current-Phase Scope

This phase is an open-data synthesis for prioritizing endothelial transition hypotheses. It does not estimate patient fate or establish that SCA changes a transition.

## Functional Requirements

- **[REQ-001]** The current-phase workflow shall construct an analysis-ready reference of adult human cardiac and pulmonary endothelial states from open processed public single-cell or single-nucleus data.
- **[REQ-002]** The current-phase workflow shall define candidate transitions as reproducible within-lineage molecular-state adjacencies.
- **[REQ-003]** WHERE eligible public directionality or perturbation data are included, the current-phase workflow shall attach that evidence only to compatible candidate adjacencies and preserve its qualification.
- **[REQ-004]** The current-phase workflow shall map public SCA-relevant endothelial and immune programs to reference states and eligible candidate adjacencies without representing a disease effect.
- **[REQ-005]** The current-phase workflow shall produce uncertainty-ranked transition hypotheses and a plan naming the minimum data needed for direct follow-up.
- **[REQ-006]** The current-phase workflow shall exclude direct patient-tissue SCA fate estimation, causal disease-shift claims, matched-SCA-tissue testing, longitudinal proof, lineage tracing, generalized platform work, controlled-access data, ligand–receptor networks, and clinical prediction.

## Acceptance Criteria

- **[AC-001]** (References REQ-001)
  - **Given** the admitted-source manifest
  - **When** the reference build is reviewed
  - **Then** every input is open and processed, and the output contains eligible adult human heart and lung endothelial states

- **[AC-002]** (References REQ-002)
  - **Given** the candidate-transition catalog
  - **When** each retained pair is inspected
  - **Then** every pair has a documented lineage, adjacency rule, and reproducibility result

- **[AC-003]** (References REQ-003)
  - **Given** a candidate adjacency with attached directional or perturbational evidence
  - **When** its evidence record is inspected
  - **Then** the record identifies the compatible edge, public source, evidence type, assumptions, and limits

- **[AC-004]** (References REQ-004)
  - **Given** the SCA-relevant mapping ledger
  - **When** a mapped program is reviewed
  - **Then** it links to a reference state or eligible adjacency and is labeled as compatibility rather than a disease effect

- **[AC-005]** (References REQ-005)
  - **Given** the final hypothesis set
  - **When** the deliverables are reviewed
  - **Then** every hypothesis has an uncertainty-qualified rank and a specific next-data requirement

- **[AC-006]** (References REQ-006)
  - **Given** the final methods, results, and claims
  - **When** the scope checklist is applied
  - **Then** none of the listed excluded activities or claims is present

# Feature 2: Data Access and Provenance

Every result needs an auditable path to an eligible public source. Missing access, terms, integrity, or lineage stops the dependent analysis.

## Functional Requirements

- **[REQ-007]** WHEN a source is proposed for admission, the source registry shall record its accession, release, retrieval URL and date, expected size, format, checksum method, terms, citation, donor and sample counts, assay, anatomy, covariates, preprocessing, and known limitations.
- **[REQ-008]** WHEN source retrieval begins, the source registry shall verify access, expected files, sizes, checksums, metadata fields, and terms before analysis.
- **[REQ-009]** IF source access, terms, retrieval, integrity, or required metadata cannot be verified, THEN the current-phase workflow shall quarantine the source and its dependent artifacts until a predeclared open substitute passes admission or the affected analysis is stopped and recorded as a gap.
- **[REQ-010]** The provenance layer shall assign stable identifiers to every source and derived artifact and record input checksums, parent artifacts, code version, parameters, label mappings, exclusions, and eligibility decisions.
- **[REQ-011]** WHEN an artifact becomes a release candidate, the provenance layer shall verify a complete trace from the artifact to its exact sources and transformations.
- **[REQ-012]** IF a release candidate lacks a complete provenance trace, THEN the provenance layer shall invalidate that artifact and its descendants until the lineage is repaired and the artifacts are rebuilt.

## Acceptance Criteria

- **[AC-007]** (References REQ-007)
  - **Given** a proposed source
  - **When** its admission record is validated
  - **Then** all specified registry fields are populated or the source is marked ineligible

- **[AC-008]** (References REQ-008)
  - **Given** a retrieval candidate and its admission record
  - **When** the retrieval preflight completes
  - **Then** the recorded access, files, sizes, checksums, metadata fields, and terms all pass before analysis starts

- **[AC-009]** (References REQ-009)
  - **Given** a source with a failed admission or retrieval check
  - **When** dependency status is inspected
  - **Then** the source and every dependent artifact are quarantined, replaced by an admitted substitute, or recorded as a stopped analysis gap

- **[AC-010]** (References REQ-010)
  - **Given** any source or derived artifact identifier
  - **When** its provenance record is queried
  - **Then** the record contains every applicable lineage field and resolves each parent identifier

- **[AC-011]** (References REQ-011)
  - **Given** a release candidate
  - **When** a result-to-source trace is executed
  - **Then** the trace resolves every parent, transformation, mapping, exclusion, and eligibility decision without a missing link

- **[AC-012]** (References REQ-012)
  - **Given** a release candidate with an incomplete trace
  - **When** the provenance control is triggered
  - **Then** the candidate and its descendants remain invalid until a rebuilt trace passes AC-011

# Feature 3: Healthy Endothelial Reference

Harmonization keeps organ biology and source context visible. Integration cannot erase distinctions that the source data support.

## Functional Requirements

- **[REQ-013]** The reference workflow shall represent eligible cardiac arterial, venous, capillary, lymphatic, and endocardial states and eligible pulmonary arterial, venous, general capillary, aerocyte, and lymphatic states.
- **[REQ-014]** WHILE constructing the healthy reference, the harmonization workflow shall preserve source labels, harmonized labels, stable gene identifiers, anatomy, donor, sample, assay, study, condition, and evidence role, use explicit gene and state crosswalks, and remove known duplicate observations.
- **[REQ-015]** WHEN a reference build completes, the harmonization workflow shall measure unresolved mappings, marker retention, donor and study leverage, batch separation, and leave-one-study-out state stability.
- **[REQ-016]** IF harmonization confounds biology with study or assay effects, erases a supported organ-specific state, or leaves an ambiguous state mapping, THEN the harmonization workflow shall split or relabel the affected state, analyze sources separately, downgrade the cross-source claim, or abstain.

## Acceptance Criteria

- **[AC-013]** (References REQ-013)
  - **Given** the admitted source coverage
  - **When** reference-state coverage is reported
  - **Then** every listed state with adequate annotation and donor support is represented and every omission has a recorded reason

- **[AC-014]** (References REQ-014)
  - **Given** a harmonized observation
  - **When** its retained metadata and mappings are inspected
  - **Then** every available specified field resolves to the source record and no known duplicate observation remains

- **[AC-015]** (References REQ-015)
  - **Given** a completed reference build
  - **When** the harmonization report is generated
  - **Then** the report contains a result for each specified diagnostic at state, donor, and study level where applicable

- **[AC-016]** (References REQ-016)
  - **Given** a failed harmonization diagnostic
  - **When** the affected state is reviewed
  - **Then** its record contains one selected response and no unsupported cross-source state claim remains

# Feature 4: Transition and SCA-Evidence Inference

Adjacency, direction, perturbation, and SCA-signal compatibility are distinct evidence fields. Cells do not substitute for donors, and association does not substitute for causation.

## Functional Requirements

- **[REQ-017]** The inference workflow shall retain a candidate state pair only when it remains within a defensible lineage, meets predeclared state-space and marker-program adjacency rules, and has the required donor or source support.
- **[REQ-018]** WHEN a candidate state pair first passes eligibility, the inference workflow shall label it as an undirected adjacency and keep adjacency, directionality, perturbation, and SCA-signal compatibility in separate fields.
- **[REQ-019]** WHERE eligible directional evidence exists, the inference workflow shall derive direction only from a compatible public velocity object, time-resolved series, or direct perturbation and record the method, source, assumptions, donor coverage, transport limits, and contradictory evidence.
- **[REQ-020]** IF a directional result cannot be reproduced or depends on ineligible or circular evidence, THEN the inference workflow shall remove the directional label and retain only the independently supported undirected adjacency.
- **[REQ-021]** WHEN an SCA-relevant program is mapped, the evidence ledger shall identify it as endothelial or immune, preserve its in-vitro, bulk, genotype, comparator, and disease-only status, and restrict its interpretation to that evidence role.
- **[REQ-022]** IF inference review detects cell-level pseudoreplication, discovery–validation reuse, genotype pooling, unsupported temporal or causal language, or interpretation beyond the source transport grade, THEN the inference workflow shall recompute at donor level, remove the arrow, downgrade the evidence, narrow the interpretation, or drop the hypothesis.

## Acceptance Criteria

- **[AC-017]** (References REQ-017)
  - **Given** a retained candidate pair
  - **When** its eligibility record is evaluated
  - **Then** the lineage, adjacency, marker-program, and donor-or-source gates all pass under their predeclared thresholds

- **[AC-018]** (References REQ-018)
  - **Given** a newly eligible candidate pair
  - **When** its record is created
  - **Then** its initial direction is undirected and the four evidence fields are independently populated or explicitly absent

- **[AC-019]** (References REQ-019)
  - **Given** a directionally labeled edge
  - **When** its directional evidence record is validated
  - **Then** the record identifies one eligible evidence class and contains every specified qualification field

- **[AC-020]** (References REQ-020)
  - **Given** a directional result that fails reproducibility or evidence eligibility
  - **When** the directionality control is applied
  - **Then** the released edge has no direction and retains only independently passing adjacency evidence

- **[AC-021]** (References REQ-021)
  - **Given** a mapped SCA-relevant program
  - **When** its evidence role is inspected
  - **Then** every applicable source-status field is explicit and the permitted interpretation matches that role

- **[AC-022]** (References REQ-022)
  - **Given** an inference record with a listed failure
  - **When** the biological-inference control is triggered
  - **Then** the record documents one permitted response and no failed claim remains in the ranked output

# Feature 5: Compute and Validation Controls

The workflow is designed for bounded, zero-budget execution on processed endothelial subsets. Scientific eligibility cannot be relaxed to solve a resource problem.

## Functional Requirements

- **[REQ-023]** WHEN a source passes admission, the compute workflow shall record storage and memory estimates, use sparse or backed endothelial-subset objects, pin the environment and random seeds, and define checksummed checkpoints.
- **[REQ-024]** WHEN a full-scale analysis is proposed, the compute workflow shall complete a small preflight and record peak memory, disk use, runtime, completion status, and repeatability.
- **[REQ-025]** IF full-scale execution exceeds available memory, disk, or session limits, THEN the compute workflow shall resume from a verified checkpoint or apply a documented biologically neutral partition and repeat scope and eligibility review before changing any cohort or endpoint.
- **[REQ-026]** The validation workflow shall reserve donor or study holdouts before discovery, bootstrap donors or studies rather than cells, and test donor leave-out, study leave-out, alternative adjacency specifications, null behavior, and directional concordance.
- **[REQ-027]** IF an edge lacks independent support, fails a stability criterion, or cannot be distinguished from study or assay effects, THEN the validation workflow shall fail its eligibility gate, label it exploratory or unresolved, and convert the missing validation into a specific acquisition requirement.

## Acceptance Criteria

- **[AC-023]** (References REQ-023)
  - **Given** an admitted source
  - **When** its execution plan is inspected
  - **Then** the plan contains both resource estimates, an endothelial-subset storage method, pinned reproducibility settings, and checksummed checkpoint locations

- **[AC-024]** (References REQ-024)
  - **Given** a proposed full-scale analysis
  - **When** preflight evidence is reviewed
  - **Then** each specified resource and repeatability measure is recorded before full-scale execution starts

- **[AC-025]** (References REQ-025)
  - **Given** an execution that exceeds an available resource limit
  - **When** processing resumes
  - **Then** it uses a verified checkpoint or documented neutral partition, and any cohort or endpoint change has a completed scope and eligibility review

- **[AC-026]** (References REQ-026)
  - **Given** the finalized validation design
  - **When** its partitions and results are audited
  - **Then** holdouts predate discovery, biological resampling uses donors or studies, and every specified sensitivity test has a recorded result

- **[AC-027]** (References REQ-027)
  - **Given** an edge that fails an independence, stability, or confounding check
  - **When** validation status is assigned
  - **Then** the edge is absent from the eligible set, marked exploratory or unresolved, and linked to a concrete missing-data requirement

# Feature 6: Ranking, Uncertainty, and Delivery

The deliverable is an inspectable evidence ranking, not a hidden score or clinical product. Failed controls and null findings remain visible.

## Functional Requirements

- **[REQ-028]** The ranking workflow shall report separate component scores for healthy-reference reproducibility, independent directional support, SCA-signal compatibility, donor and study coverage, analysis sensitivity, and biological transport distance.
- **[REQ-029]** The ranking workflow shall report adjacency confidence, directionality uncertainty, mapping uncertainty, transport limitations, counterevidence, permitted interpretation, and triggered-control history for every hypothesis.
- **[REQ-030]** WHEN a release candidate is prepared, the release workflow shall verify links, manifests, checksums, exclusions, null results, repository status, independent reproduction, and causal and clinical claim boundaries.
- **[REQ-031]** IF a release check fails, THEN the release workflow shall withhold the affected artifact, correct and version it, record any changed rank, and retain the failed result in the decision log.
- **[REQ-032]** The release workflow shall publish a source registry and provenance manifest, the harmonized reference, the adjacency catalog, qualified directional evidence, the SCA mapping ledger, ranked hypothesis records, a failure and decision log, and the prioritized follow-up data plan.
- **[REQ-033]** IF any prevention or detection control fails, THEN the current-phase workflow shall record an exclusion, downgrade, abstention, correction, stopped analysis, or explicit follow-up requirement instead of an undocumented workaround.
- **[REQ-034]** The release workflow shall describe group-level molecular evidence, distinguish genotype, ancestry, geography, environment, treatment access, and self-identified race, disclose representation gaps, and avoid individual prediction or disease-reductive language.

## Acceptance Criteria

- **[AC-028]** (References REQ-028)
  - **Given** a ranked hypothesis
  - **When** its ranking record is inspected
  - **Then** all six component scores are visible with their source evidence and no unexplained aggregate rank is the sole result

- **[AC-029]** (References REQ-029)
  - **Given** a ranked hypothesis
  - **When** its uncertainty record is inspected
  - **Then** every specified uncertainty, limitation, counterevidence, interpretation, and control field is populated or explicitly marked absent

- **[AC-030]** (References REQ-030)
  - **Given** a release candidate
  - **When** the release checklist is executed
  - **Then** every specified check has a passing result before publication

- **[AC-031]** (References REQ-031)
  - **Given** a failed release check
  - **When** release status is reviewed
  - **Then** the affected artifact is unpublished until its corrected version, changed rank, and retained failure record are present

- **[AC-032]** (References REQ-032)
  - **Given** a completed current-phase release
  - **When** its deliverable manifest is validated
  - **Then** every specified artifact is present, versioned, and linked from the manifest

- **[AC-033]** (References REQ-033)
  - **Given** any failed prevention or detection control
  - **When** its control record is audited
  - **Then** exactly one documented response category is present and no undocumented substitution or workaround is used

- **[AC-034]** (References REQ-034)
  - **Given** the released narrative and result cards
  - **When** the ethics and interpretation checklist is applied
  - **Then** all specified distinctions and disclosures are present and no individual, clinical, or disease-reductive claim remains
