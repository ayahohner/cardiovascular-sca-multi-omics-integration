# Current phase: open-data endothelial transition synthesis

**Project:** Cardiovascular SCA Multi-Omics Integration · **Status:** Current implementation scope · **Version:** 0.3 · **Revised:** 2026-08-02

## Aim

Build a reproducible healthy heart-and-lung endothelial state reference, define biologically adjacent candidate transitions within that reference, and use only public evidence to rank which transitions deserve direct testing in SCA tissue.

The output is a set of **transition hypotheses**, not an estimate of what patient cells did or will do. “Transition” means a reproducible relationship between nearby molecular states. It becomes directional only when an independent public time-resolved, velocity, or perturbation dataset supports direction. Mapping an SCA-relevant signal to a transition means molecular compatibility, not proof that SCA changes its probability.

## Scope contract

### In scope now

1. **Healthy endothelial reference.** Build an analysis-ready reference of adult cardiac and pulmonary endothelial states from open, processed public single-cell or single-nucleus data.
2. **Candidate transitions.** Define reproducible, biologically adjacent state pairs using within-lineage expression structure, marker programs, anatomy, and donor or study support.
3. **Supported directionality.** Add direction or perturbation evidence only where suitable public data exist; otherwise retain an undirected adjacency.
4. **SCA-relevant signal mapping.** Map public endothelial perturbation signatures and donor-level immune programs associated with SCA onto states and candidate transitions, preserving each source's evidence role.
5. **Prioritized follow-up.** Rank hypotheses with explicit uncertainty, counterevidence, transport limitations, and a concrete data-acquisition or experimental follow-up plan.

### Out of scope or deferred

- direct patient-tissue SCA fate estimation;
- claims that SCA causes a transition or changes its probability;
- testing in matched SCA heart or lung tissue until such data are openly available;
- longitudinal proof or lineage tracing;
- generalized multi-modal data-plane or platform construction;
- controlled-access datasets;
- ligand–receptor or cell–cell communication networks; and
- diagnosis, prognosis, treatment recommendation, or other clinical prediction.

These boundaries govern implementation. Any expansion requires a separate scope decision.

## Scientific question

**Which endothelial state adjacencies are reproducible across healthy adult heart and lung references, which have credible public evidence for direction, and which state or transition programs are most compatible with public SCA-relevant endothelial and immune signals?**

This question supports hypothesis prioritization. It cannot establish a disease-induced shift because the current phase has no matched SCA cardiopulmonary tissue, longitudinal observation, or lineage measurement.

## Biological boundary

The modeled populations are adult human endothelial states from heart and lung:

- cardiac arterial, venous, capillary, lymphatic, and endocardial endothelium;
- pulmonary arterial and venous endothelium, general capillary cells, aerocytes, and lymphatic endothelium; and
- additional endothelial subtypes only when their source annotations and donor support are adequate.

Mural, stromal, cardiac, pulmonary, and immune populations may provide annotation context. Immune data enter as external SCA-relevant programs; this phase does not infer intercellular signaling edges.

## Evidence design

### 1. Construct the healthy reference

Use processed matrices and source metadata rather than reprocessing raw reads. Preserve donor, tissue, anatomical region, assay, study, original labels, and harmonized labels. Align comparable endothelial states across heart and lung without forcing organ-specific or outlying states to merge.

### 2. Define candidate transitions

A candidate pair must remain within a defensible lineage, be adjacent under predeclared state-space and marker-program criteria, and meet minimum donor or source support. Recompute adjacency under reasonable integration, feature-selection, and neighborhood choices. An edge is initially an undirected biological adjacency, not a trajectory.

### 3. Attach directionality or perturbation evidence

Direction may be added from a public velocity object with suitable spliced and unspliced information, a time-resolved series, or a direct perturbation whose biological transport is stated. Each edge records the method, source, assumptions, and contradictory evidence. If direction cannot be reproduced, the edge remains undirected.

### 4. Map SCA-relevant signals

Endothelial and immune signatures are scored against state programs and, where meaningful, the expression difference between adjacent states. Evidence roles remain separate:

- an *in vitro* endothelial perturbation is a perturbational prior, not patient tissue;
- donor-replicated PBMC data support systemic immune association, not endothelial identity; and
- a disease-only cohort without a healthy comparator supplies context, not a disease effect.

### 5. Rank hypotheses

Rank each transition by healthy-reference reproducibility, strength and independence of directionality evidence, compatibility with SCA-relevant signals, donor and study coverage, sensitivity to analysis choices, and biological transport distance. Publish the component scores and evidence grade; do not hide them inside a single unexplained rank.

## Candidate open sources

All sources require a retrieval, license, metadata, and checksum audit before analysis. Only open processed data are eligible for this phase.

| Source | Current-phase role | Evidence limit |
| --- | --- | --- |
| [Human Heart Cell Atlas v1](https://www.heartcellatlas.org/v1.html) | Healthy cardiac endothelial states and donor/anatomical context | Healthy reference, not an SCA control cohort |
| [Integrated Human Lung Cell Atlas v1.0](https://data.humancellatlas.org/hca-bio-networks/lung/atlases/lung-v1-0) | Healthy pulmonary endothelial states and cross-study context | Integrated reference with study and mapping effects |
| Heart Cell Atlas vascular velocity object, if its files and metadata pass audit | Candidate directional evidence within healthy cardiac vasculature | Direction is method-dependent and covers only eligible cells and donors |
| [GEO GSE274340](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274340) | Descriptive endothelial perturbation signature | Two-condition HUVEC model; not patient tissue and not donor-replicated disease evidence |
| [GEO GSE53441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53441) | Donor-level HbSS-versus-HbAA PBMC program | Bulk, older-platform, systemic immune evidence |
| [GEO GSE84632](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84632) | Optional SCD PBMC context | No internal healthy comparator; not clean replication of a disease effect |

The controlled-access EGA HSPC source described in the earlier proposal is excluded from this phase. So are unverified or unavailable SCA heart, lung, spatial, longitudinal, chromatin, proteomic, and clinical datasets.

## Provenance and reproducibility

Each source must have a versioned record containing its accession, retrieval URL and date, checksum, license or terms, citation, file format, donor and sample counts, assay, anatomy, available covariates, source preprocessing, exclusions, and known limitations. Raw source artifacts remain immutable; every derived object records its input checksums, code version, parameters, label mappings, and eligibility decision.

Common semantics are deliberately narrow: stable gene identifiers, source and harmonized endothelial labels, anatomy, donor, sample, assay, study, condition, and evidence role. This phase will create the manifests and analysis objects it needs, not a generalized cohort service or multi-modal platform.

Donor is the biological replication unit when donor identities exist. Cells are not treated as independent people. Sources without biological replication can contribute descriptive signatures but not population-level inferential claims.

## Staged execution

1. **Audit and freeze inputs.** Verify open access and terms; record provenance; define eligibility, adjacency, directionality, mapping, ranking, abstention, and failure-response rules before modeling.
2. **Build and validate the reference.** Harmonize eligible heart and lung endothelial states while retaining organ, donor, assay, study, and source-label structure.
3. **Generate candidate adjacencies.** Construct and stress-test undirected state pairs; discard edges that depend on one donor, study, or fragile parameter choice.
4. **Add public evidence.** Attach qualified directionality or perturbation support, then map SCA-relevant endothelial and immune programs without collapsing evidence types.
5. **Rank and report.** Produce evidence-graded hypotheses, uncertainty and sensitivity summaries, unresolved questions, and a prioritized data-acquisition or experimental plan.

## Prospective postmortem and controls

Assume the phase failed to produce trustworthy, reproducible hypotheses. The most foreseeable causes and required controls are:

| Failure mode | Prevention | Detection | Response |
| --- | --- | --- | --- |
| **Data access:** a link moves, a processed object is unavailable or impractically large, a portal throttles retrieval, or terms are unclear. | Admit a source only after recording a stable accession, release, URL, expected size, format, checksum method, terms, and a processed-data retrieval path. Cache one immutable copy and use endothelial subsets rather than raw reads or whole-atlas reprocessing. | Run a retrieval and integrity preflight before analysis; compare expected files, sizes, checksums, and metadata fields with the source record. | Quarantine the source and all dependents. Retry conservatively or use a predeclared open mirror or subset; otherwise record the gap and stop the affected analysis. Never substitute silently. |
| **Provenance:** a result cannot be traced to an exact source, transformation, label map, or eligibility decision. | Give every source and derived artifact a stable identifier; record input checksums, code version, parameters, mappings, exclusions, and parent artifacts. | Sample result-to-source traces during each stage and require a complete trace for every release candidate. | Invalidate the untraceable artifact and everything derived from it, repair the lineage record, and rebuild before release. |
| **Harmonization:** gene or state mappings are wrong, duplicate observations inflate support, study effects masquerade as biology, or integration erases organ-specific states. | Preserve source labels and donor, organ, assay, and study fields; use explicit gene and state crosswalks; deduplicate known reused observations; align only eligible endothelial strata. | Inspect unresolved mappings, source-specific marker retention, donor and study leverage, batch separation, and leave-one-study-out state stability. | Split or relabel the affected state, analyze sources separately, downgrade cross-source claims, or abstain when biology and study cannot be separated. |
| **Biological inference:** adjacency is reported as temporal fate; cells are treated as independent donors; direction is circular; or HUVEC, PBMC, genotype, and disease-only evidence are overgeneralized. | Keep adjacency, directionality, perturbation, and SCA-signal compatibility as separate fields. Use donor-level replication, predeclare evidence roles, and forbid disease-shift or causal language. | Audit for cell-level pseudoreplication, unsupported arrows, reused discovery and validation data, genotype pooling, and conclusions stronger than the source's transport grade. | Remove the arrow, recompute at donor level, downgrade to descriptive evidence, narrow the interpretation, or drop the hypothesis. |
| **Compute:** full objects exceed memory or disk, free sessions expire, or stochastic runs are not reproducible. | Estimate storage and memory before retrieval; use sparse, backed, and endothelial-subset objects; pin environments and random seeds; checkpoint derived artifacts with checksums. | Run a small preflight and record peak memory, disk, runtime, completion status, and repeatability before the full analysis. | Partition or reduce only through a documented, biologically neutral rule; resume from verified checkpoints. If feasibility changes the cohort or endpoint, rerun scope and eligibility review. |
| **Validation:** too few independent donors or studies remain, holdouts are contaminated, or transitions disappear under reasonable specifications. | Reserve donor or study holdouts before discovery, keep validation sources independent where possible, and predeclare stability and abstention criteria. | Perform donor and study leave-out analyses, alternative adjacency specifications, null checks, and directionality concordance tests. | Fail the edge's eligibility gate, label it exploratory or unresolved, and convert the missing validation into a specific acquisition requirement. |
| **Delivery:** ranks are opaque, negative results vanish, links break, raw data enter Git, or wording implies clinical or causal validity. | Use a fixed result card with component evidence, uncertainty, counterevidence, provenance, permitted claim, and next-data requirement. Keep raw data outside Git and run a release checklist. | Require independent reproduction and claim-boundary review; verify links, manifests, checksums, exclusions, null results, and repository status. | Withhold the release, correct and version the affected artifact, publish the correction and changed rank, and retain the failed result in the decision log. |

A control failure is an analysis result. It must lead to exclusion, downgrade, abstention, or an explicit follow-up requirement—not an undocumented workaround.

## Validation and uncertainty

- Hold out donors and, where possible, an entire source study when evaluating reference stability.
- Rebuild state adjacency under alternative integration, feature, and neighborhood specifications.
- Bootstrap donors or studies rather than cells for biological uncertainty.
- Record mapping uncertainty for ambiguous or organ-specific state labels.
- Separate directionality uncertainty from adjacency confidence.
- Treat transport from HUVEC or PBMC evidence to adult cardiac or pulmonary endothelium as a named limitation.
- Abstain when an edge lacks adequate donor support, depends on one analysis choice, or cannot be distinguished from study or assay effects.

No interval or score should imply observed temporal fate. A null result—no transition meeting the evidence threshold—is valid and should sharpen the follow-up data plan.

## Outputs

1. A versioned open-source registry, provenance manifest, and source-admission record.
2. A harmonized healthy heart-and-lung endothelial state reference with source labels retained.
3. A reproducible catalog of candidate state adjacencies and any separately qualified directional evidence.
4. A mapping ledger for public SCA-relevant endothelial and immune programs.
5. Ranked transition hypotheses with component evidence grades, uncertainty, counterevidence, permitted interpretation, and triggered-control history.
6. A failure and decision log plus a prioritized plan for the matched tissue, longitudinal, lineage, perturbation, or validation data needed next.

## Ethics and interpretation

The phase uses public research data and produces group-level hypotheses, not individual predictions. Genotype, ancestry, geography, environment, treatment access, and self-identified race remain distinct concepts; race is not used as a biological proxy. Dataset representation and missing context must appear beside results. Language should describe evidence about molecular programs without reducing people with SCA to a disease label.

## Success criteria

The phase succeeds when:

- another analyst can recreate the eligible inputs and reference from the recorded accessions, checksums, mappings, and parameters;
- retained adjacencies are stable across donors or studies and reasonable analysis choices;
- directional labels cite independent eligible evidence and otherwise remain absent;
- SCA-relevant mappings preserve their endothelial, immune, *in vitro*, bulk, and disease-only evidence roles;
- every ranked hypothesis includes uncertainty, transport limitations, counterevidence, and a bounded interpretation;
- every triggered control has a recorded exclusion, downgrade, abstention, correction, or follow-up response; and
- the final acquisition plan names the minimum data needed to test the leading hypotheses directly.
