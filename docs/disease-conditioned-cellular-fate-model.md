# A unified data plane for disease-conditioned cellular fate across the cardiopulmonary vascular axis

**Status:** Sample project specification — planning only · **Version:** 0.1 · **Prepared:** 2026-08-01

## The pitch

Sickle cell anemia begins with altered hemoglobin, but its burden is expressed through a living vascular system: red cells, endothelium, mural cells, immune cells, stroma, lung, and heart continually act upon one another. Healthy cell atlases give us increasingly precise maps of those tissues. Sickle-cell datasets, by contrast, capture scattered views of endothelial stress, circulating immunity, and hematopoietic programming.

The immediate bottleneck is not finding more files. It is making unlike evidence comparable without stripping away the facts that determine what it means: anatomical site, cell semantics, assay behavior, donor and cohort context, processing history, quality, consent, and access conditions. A model built directly on loosely aggregated studies could learn the difference between laboratories and call it disease.

The project's primary product is therefore a **unified, analysis-ready data plane** for healthy-versus-SCA cardiovascular cell-state research. “Unified” describes a common semantic and analytical interface, not one physical database and not permission to copy controlled data into an open store. The disease-conditioned fate model is the first demanding consumer and proof of value: it will estimate how the *probabilities* of candidate cell-state transitions and communication relationships differ between a healthy reference and SCA. It will not claim to replay a patient's disease, prove a mechanism, or infer literal lineage from cross-sectional samples. Its value is a transparent, uncertainty-aware map of where disease-associated evidence most strongly reweights the healthy cellular landscape—and an enduring data foundation on which later questions can be asked without rebuilding context each time.

## Precise research question

Among defined endothelial, mural, stromal, cardiac, pulmonary, and immune populations along the adult human cardiopulmonary vascular axis, **which candidate cell-state transitions and intercellular signaling relationships have different estimated probabilities in SCA-associated data than in harmonized healthy references, after accounting for donor, tissue, assay, and available clinical covariates?**

Here, “fate” means a modeled movement between molecularly defined cell states. A transition such as quiescent-to-activated endothelium is a probabilistic hypothesis supported by state similarity, directionality evidence when available, and external biology. It is not lineage tracing. A disease-associated probability shift is an association conditional on the included data and model; it is not causal proof that SCA produced the shift.

The primary disease group is SCA as defined by each source. HbSS should be the v1 anchor when genotype is available. HbSβ⁰ and other sickling genotypes must be retained as distinct strata rather than silently pooled.

## Primary value proposition: an analysis-ready data plane

The data plane is a governed analytical abstraction across distributed sources. It does more than inventory datasets. For every cell, sample, donor, cohort, tissue region, measurement, and derived feature that can legally be exposed, it supplies a stable meaning, context, provenance chain, quality status, and intended-use boundary.

An analysis request should return an **evidence bundle**, not a naked matrix. The bundle contains the modality-appropriate data view; donor and cohort context; original and harmonized labels; anatomy and condition terms; source and transformation lineage; quality tier and failed checks; access and consent restrictions; missingness; and donor-grouped analysis splits. A researcher should be able to ask “what supports this transition?” and walk from a result back to the permitted source records and every transformation between them.

| Data-plane contract | What it guarantees |
| --- | --- |
| Common semantics | Stable identifiers and explicit crosswalks for genes, cell types and states, anatomy, disease genotype, phenotype, assay, donor, and cohort—without discarding source labels. |
| Modality awareness | scRNA-seq, snRNA-seq, spatial assays, bulk expression, proteins, imaging, functional measurements, and clinical variables are linked where defensible but are not forced into one false-equivalent matrix. |
| Cohort and context | Donor, treatment, transfusion, clinical state, age, sex as recorded, genotype, geography, collection site, time point, tissue, and batch travel with every analytical view. |
| End-to-end provenance | Accession, version, consent and license terms, checksums, retrieval, preprocessing, harmonization, exclusions, parameters, and code version remain queryable. |
| Quality and eligibility | Each source or subset has a visible quality tier, analysis-specific eligibility flags, and reasons for exclusion or abstention. |
| Governed interoperability | Open and controlled-access data can participate in one logical analysis plan while physical data, row-level metadata, and exported derivatives remain inside their permitted boundaries. |
| Analysis interface | Downstream models receive repeatable cohort definitions, comparable features, donor-aware splits, spatial relations, and provenance—not ad hoc file paths. |

The data plane is successful when two analysts can express the same cohort and evidence question, receive equivalent permitted inputs, and understand every remaining mismatch. It is not successful merely because files share a folder or object format.

## Deliberate v1 scope

### Biological boundary

V1 covers the **cardiopulmonary vascular axis**: cardiac chambers and myocardium, coronary and intramyocardial vessels, pulmonary arteries and microvasculature, lung vascular niches, and the circulating blood and immune compartments that interact with them. This boundary is broad enough to connect oxygenation, vascular stress, inflammation, and cardiac adaptation, yet narrow enough to resist becoming a whole-body SCA atlas.

The default comparison is healthy adult biology versus clinically stable SCA. Acute pain crisis, acute chest syndrome, pulmonary hypertension, transfusion state, hydroxyurea or hydroxycarbamide exposure, age, and other clinical contexts are covariates or separately labeled strata—not interchangeable examples of one “sick” condition. Pediatric and adult samples should not be pooled without a predeclared age model.

### Cell populations in scope

| Compartment | Populations | Place in the model |
| --- | --- | --- |
| Cardiac vascular | Arterial, venous, capillary, lymphatic, and endocardial endothelial states | Primary state-transition targets |
| Pulmonary vascular | Arterial and venous endothelial cells, general capillary cells, aerocytes, and lymphatic endothelial cells | Primary state-transition targets |
| Mural and stromal | Pericytes, vascular smooth-muscle cells, fibroblasts, myofibroblast-like states, and mesothelial cells where relevant | Vessel stability, remodeling, and signaling context |
| Immune | Classical and non-classical monocytes, macrophages, neutrophils, dendritic cells, T cells, B cells, and natural killer cells | Circulating and tissue inflammatory programs |
| Blood and hematopoietic | Erythrocytes, reticulocytes, platelets, and hematopoietic stem/progenitor cells | Perturbation inputs and upstream context; mature red cells are not expected to be well represented by scRNA-seq |
| Organ context | Cardiomyocytes, alveolar epithelial cells, and other abundant resident populations | Receivers, senders, and spatial anchors rather than exhaustive fate targets in v1 |

Rare cell types remain visible, but no transition claim should be made for a population without adequate donor support.

## Conceptual architecture

The architecture separates source truth, harmonization, eligibility, access, and modeling so that biological and governance assumptions do not disappear inside one embedding.

1. **Source registry and policy envelope.** Register every open or controlled source before retrieval. Capture stable identifiers, versions, cohort definitions, consent and data-use terms, custodians, permitted purposes, and permitted derivatives. A source can be discoverable without being accessible.
2. **Modality-specific landing zones.** Preserve immutable source files and source-native metadata. scRNA-seq, snRNA-seq, spatial, bulk, protein, imaging, functional, and clinical inputs receive modality-appropriate validation rather than a universal cleaning recipe.
3. **Common semantic plane.** Resolve genes, cell types and states, anatomy, genotype, phenotype, assay, sample, donor, and cohort to versioned identifiers. Retain original labels beside harmonized terms and record every exact, broader, narrower, or unresolved mapping.
4. **Modality-aware harmonization.** Normalize and integrate within defensible comparison sets, preserve raw measurements, quantify batch and study effects, and connect modalities through declared mappings or shared features. Healthy and SCA evidence are aligned to a reference coordinate system without erasing condition-specific or out-of-reference states.
5. **Cohort and context plane.** Attach clinical state, treatment, transfusion, age, sex as recorded, genotype, geography, anatomical site, time point, collection protocol, and technical batch at their correct donor, sample, or measurement level. Missing, withheld, inapplicable, and unknown values remain distinct.
6. **Quality and eligibility plane.** Assign quality tiers and analysis-specific eligibility at source, donor, sample, cell, and feature level. Quality is not a single score: a cohort may be suitable for program replication yet unsuitable for transition inference.
7. **Governed analysis interface.** Materialize consent-aware evidence bundles for cohort selection, reference mapping, donor-level comparisons, spatial neighborhoods, state graphs, and validation splits. The interface exposes provenance and policy with the data and can execute within a controlled environment when records may not leave it.
8. **Disease-conditioned model.** The healthy atlas view supplies baseline state and communication graphs. SCA evidence estimates perturbation terms that can increase, decrease, or leave unchanged transition probabilities and signaling edges, conditional on measured context. The model writes results back with the exact evidence bundle, quality tier, and model version that produced them.

In compact form, the fate model asks for the probability of a candidate next state given a current state, anatomical context, donor context, and condition. Disease is a perturbation of the healthy transition system—not a separate universe and not a label the model is trained merely to classify. The data plane makes each conditioning variable and evidence boundary explicit enough for that comparison to be audited.

No single trajectory method is privileged in the specification. RNA velocity is eligible only for datasets with suitable unspliced counts and defensible kinetics. Pseudotime, optimal transport, or graph-based directionality may support a transition, but cross-sectional agreement among methods still does not turn association into observed temporal fate.

## Data plan

### Required modalities

- **Single-cell or single-nucleus RNA sequencing:** the core vocabulary of cell identity and state in healthy heart, lung, vasculature, blood, and SCA-associated samples.
- **Spatial transcriptomics or imaging-based spatial assays:** anatomical localization, neighborhood composition, and constraints on ligand–receptor claims. Healthy spatial references are useful but cannot substitute for SCA tissue.
- **Chromatin accessibility or joint multi-omics, when available:** regulatory support for state distinctions and transition direction. This is desirable, not required for the MVP.
- **Bulk transcriptomics:** donor-level immune or endothelial modules from larger SCA cohorts. Bulk data should validate or constrain cell-state programs, not be presented as cell-resolved evidence.
- **Protein and functional measurements:** cytokines, adhesion markers, flow cytometry, red-cell deformability, hemolysis-related measurements, oxygen tension, and vascular phenotypes where available. These provide orthogonal evidence and perturbation context.
- **Clinical and technical metadata:** genotype, age, sex as recorded, treatment and transfusion exposure, clinical state at sampling, cardiopulmonary phenotype, tissue and anatomical site, collection protocol, library chemistry, sequencing depth, and batch.

Harmonization must respect what each modality measures. Single-cell and single-nucleus assays require distinct capture and cell-composition checks; spatial spots or segmented cells carry coordinates and neighborhood uncertainty; bulk assays support donor-level programs but cannot be converted into invented cells; protein and functional measurements retain units, detection limits, and specimen context. Cross-modality links should be typed as direct measurement, aggregation, mapping, imputation, or prior so downstream analyses can include or reject them deliberately.

### Candidate public-data snapshot

The following are **candidates, not a locked bill of materials**. Availability, consent, licensing, cohort definitions, and file integrity must be re-audited before analysis. Portal facts below were checked on 2026-08-01.

| Candidate source | Proposed role | Qualification |
| --- | --- | --- |
| [Human Heart Cell Atlas v1](https://www.heartcellatlas.org/v1.html) | Healthy cardiac and vascular reference | The portal reports 486,134 cells and nuclei from six cardiac regions and 14 adult donors. It is a reference cohort, not a matched control for SCA. |
| [Integrated Human Lung Cell Atlas v1.0](https://data.humancellatlas.org/hca-bio-networks/lung/atlases/lung-v1-0) | Healthy pulmonary and endothelial reference | The portal describes an integrated atlas of more than two million cells; its healthy core spans 107 individuals. Integration and donor composition must be carried into uncertainty analyses. |
| [GEO GSE274340](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274340) | Endothelial perturbation prior | Public scRNA-seq of HUVECs exposed to healthy red cells or a mixture containing pharmacologically dehydrated red cells in a microfluidic device. It is a two-condition *in vitro mimic*, not patient endothelium or a general SCA tissue atlas. |
| [GEO GSE53441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53441) | Donor-level immune-expression comparison | Public PBMC microarray data from 24 steady-state HbSS participants and 9 healthy HbAA participants. It offers a direct comparison but is bulk, older-platform, and not cardiopulmonary tissue. |
| [GEO GSE84632](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84632) | Larger SCD PBMC replication cohort | Public PBMC microarray data from 172 adults with SCD. It can test donor-level modules, but it is not single-cell and does not provide its own healthy comparator. |
| [EGA EGAS50000001046](https://www.ega-archive.org/studies/EGAS50000001046) | Exploratory SCA hematopoietic-state evidence | The catalogue describes peripheral-blood CD34+ scRNA-seq from four people with SCA and three healthy donors. The data are controlled-access, small, and upstream of the vascular tissues of interest. |

This snapshot reveals the central gap: readily discoverable SCA data do not form a matched, spatially resolved heart–lung–blood cohort. V1 must therefore distinguish direct disease evidence, in vitro perturbational evidence, healthy spatial constraints, and transferred priors in every result. A deeper accession audit may find additional eligible sources; absence from this short list is not evidence that a dataset does not exist.

### Common semantics and provenance

Every source receives a versioned dataset card containing its accession and DOI, retrieval date, checksum, license or data-use terms, consent and access class, cohort definition, donor counts, assay, genome build, preprocessing history, clinical variables, exclusions, and known limitations. Raw files remain immutable. Every derived object records the source objects, code version, parameters, semantic mappings, quality decisions, and transformation chain that produced it.

The shared semantic model distinguishes person, donor, specimen, sample, library, assay run, cell or nucleus, spatial observation, feature, cohort, and derived result. It also separates diagnosis from genotype, clinical state from long-term disease label, treatment exposure from treatment response, and biological absence from missing or withheld data. Those distinctions prevent a convenient table join from becoming a false scientific claim.

Genes should map through stable Ensembl identifiers with HGNC symbols retained for display. Cell types, anatomy, disease, and phenotypes should use versioned terms from established ontologies where practical, with source labels preserved beside harmonized labels. Each crosswalk records whether a mapping is exact, broader, narrower, composite, or unresolved. No label is silently merged.

### Interoperability and analysis interface

Interchange should use broadly readable objects such as H5AD for cell-by-feature matrices and Parquet or TSV for metadata, but file format is only the transport layer. The durable contract is the relationship among those objects: stable entity identifiers, units, schemas, ontology versions, spatial coordinate systems, feature mappings, and provenance links.

The analysis interface should support five practical operations without prescribing a software product:

1. **Cohort builder:** select donors and samples by tissue, condition, genotype, clinical context, modality, access class, and quality eligibility while exposing overlap and confounding.
2. **Comparable feature view:** retrieve source-native measurements and harmonized features with modality, transformation, and missingness flags.
3. **Reference and spatial view:** map cells or programs to healthy reference states and retrieve anatomical or neighborhood relations with mapping confidence.
4. **Model-ready view:** emit donor-grouped matrices, covariates, allowed transition graphs, validation splits, and abstention flags for the fate model.
5. **Evidence trace:** walk from any estimate or graph edge back through transformations, quality checks, semantic crosswalks, source records, and permitted source metadata.

The same query specification should be replayable against an updated data-plane version, producing an explicit diff in included evidence rather than a silently changed cohort.

### Quality tiers

Quality tiers are cumulative gates, assigned to each relevant subset rather than stamped once on an entire study.

| Tier | Meaning | Permitted use |
| --- | --- | --- |
| Q0 — Catalogued | Identity, ownership, access class, and basic metadata are known. | Discovery and gap analysis only. |
| Q1 — Acquired and permitted | Files, checksums, consent or license terms, and local authorization are verified. | Source-level audit inside the permitted boundary. |
| Q2 — Modality-valid | Modality-specific integrity and quality checks pass, with exclusions recorded. | Within-source summaries and quality reporting. |
| Q3 — Semantically harmonized | Required identifiers, context, units, and label crosswalks meet the common contract. | Cross-source descriptive comparison with visible mapping uncertainty. |
| Q4 — Analysis-eligible | Confounding, donor support, directionality requirements, and endpoint-specific checks pass. | Named downstream analysis only; eligibility for one analysis does not imply eligibility for all. |

An excluded source remains catalogued with its failed gates and reason. No automated repair may promote a dataset without recording the change and rerunning the affected checks.

### Controlled-access and consent boundaries

Donor is the unit of biological replication. Cells are observations within donors, not independent people. Open and controlled-access data remain physically separated according to their terms. The logical data plane may dispatch an approved analysis into a secure environment and return only permitted aggregates; it must not use harmonization as a route around consent, data-use limits, export review, or minimum-cell and minimum-donor disclosure rules.

Each evidence bundle carries purpose restrictions, export class, expiry or review dates where applicable, and the policy decision that allowed its use. Revocation or changed terms must invalidate affected derived views. Only approved, non-identifying derivatives may cross an access boundary, and a public result must remain interpretable even when its controlled inputs cannot be redistributed.

### Likely data gaps

- Little or no matched SCA single-cell and spatial profiling of human heart and lung tissue in the initial candidate set.
- Sparse longitudinal sampling, making temporal transition estimates largely model-based.
- Uneven capture of granulocytes, platelets, erythrocytes, and fragile endothelial states across common assays.
- Limited donor counts and incomplete genotype, treatment, transfusion, crisis-state, and cardiopulmonary phenotype metadata.
- Assay, tissue-processing, geography, age, sex, and ancestry distributions that differ between healthy atlases and SCA cohorts.
- Weak correspondence between HUVEC perturbation models and adult cardiac or pulmonary microvascular endothelium.
- Controlled-access delays and licenses that may prevent a uniform fully open release.

These gaps are feasibility criteria, not inconveniences to be corrected away. If the relevant condition and tissue are confounded with study or assay, the model must report non-identifiability rather than assign the difference to disease.

## Staged analysis plan

### Stage 0 — Contracts, feasibility, and preregistration

Define the semantic model, provenance events, context fields, quality gates, access-policy fields, and evidence-bundle interface before moving data. Build the candidate registry, data-use review, cohort coverage matrix, and confounding map. Freeze the primary comparison, exclusions, covariates, endpoints, and uncertainty thresholds. Stop or narrow the project if condition is inseparable from assay, tissue, or donor group.

### Stage 1 — Register and govern sources

Verify source identity, versions, checksums, consent and access conditions, permitted derivatives, and custodians. Land open and controlled files only in their approved environments. Produce source-native cohort and modality profiles without cross-study integration.

### Stage 2 — Modality-aware quality control and harmonization

Run modality-specific integrity and quality checks, then map identifiers, units, cell labels, anatomy, genotype, phenotype, and cohort context to the common semantics. Preserve source labels and unresolved mappings. Promote subsets through Q0–Q3 only when their recorded gates pass.

### Stage 3 — Materialize the analysis interface and healthy reference

Implement the logical cohort, feature, spatial, model-ready, and evidence-trace views against the harmonized sources. Assemble the healthy heart and lung coordinate system, quantify healthy donor variation, and establish spatial neighborhoods. Hold out entire donors and, where possible, an entire source study for validation.

### Stage 4 — Map disease evidence and determine eligibility

Map each SCA source independently to the frozen healthy reference. Measure out-of-reference behavior and permit an “unresolved” state. Estimate donor-level state abundance and program scores, then assign analysis-specific Q4 eligibility. Bulk, in vitro, blood, and controlled-access evidence retain distinct evidence roles.

### Stage 5 — Fit candidate state and communication graphs

Build lineage-specific graphs only where state density and donor support are sufficient. Compare graph constructions and directionality assumptions. Estimate healthy baseline probabilities and SCA-associated perturbations with hierarchical donor effects and measured covariates. Build sender–receiver hypotheses only with expression support in both populations and spatial plausibility where available.

### Stage 6 — Validate, package, and release

Run internal, external, and biological validation; test replayability of the evidence bundles; complete dataset, interface, and model cards; and freeze figures and tables from reproducible outputs. Release only derivatives allowed by source terms. Document negative and indeterminate results alongside positive findings.

## Validation plan

### Data-plane validation

- Validate schemas, stable identifiers, referential integrity, units, ontology versions, and required context at ingestion and after every transformation.
- Recompute a fixed set of cohort queries after source or contract updates and explain every membership change.
- Trace sampled model inputs and outputs end to end, including secure-environment and export-policy decisions.
- Test that withheld, unknown, absent, and not-applicable values remain distinguishable through the analysis interface.
- Confirm that open queries cannot expose controlled row-level records or disallowed small aggregates.
- Have a second analyst reproduce the primary evidence bundle from its declarative cohort and eligibility specification.

### Internal validation

- Hold out donors—not cells—for cross-validation.
- Test label transfer and transition estimates on a held-out atlas or source study.
- Bootstrap donors and repeat the analysis under alternative integration, state-definition, and trajectory choices.
- Use condition-label permutations at donor level as a null test.
- Run leave-one-donor and leave-one-study-out influence analyses.
- Check whether treatment, age, sex, genotype, assay, or sampling state explains an apparent disease shift.

### External and orthogonal validation

- Require prespecified endothelial or immune modules to show compatible direction in at least one independent SCA cohort.
- Compare cell-resolved programs with bulk PBMC, protein, flow-cytometry, red-cell functional, or clinical measurements when available.
- Test whether proposed communication edges are anatomically plausible in independent spatial or imaging references.
- Ask cardiovascular, pulmonary vascular, hematology, and patient-domain experts to review interpretation and missing context before public claims.

Wet-lab confirmation would be the logical next project, not a hidden completion criterion for this computational v1.

## Outputs

1. A versioned data-plane specification: semantic contracts, provenance model, quality tiers, policy fields, and analysis-interface contract.
2. A source registry and cardiopulmonary coverage/gap map spanning open and controlled-access candidates.
3. Replayable, consent-aware evidence bundles for the primary healthy and SCA cohorts.
4. A harmonized healthy cardiopulmonary vascular reference with traceable original labels and modality context.
5. Donor-aware maps of SCA evidence onto that reference, including unresolved regions.
6. Transition-difference tables and visual state graphs with confidence or credible intervals.
7. A spatially constrained cell-communication hypothesis network.
8. Short result cards stating evidence, counterevidence, uncertainty, quality tier, and permitted interpretation for each prioritized transition or edge.
9. A reproducible technical report and an accessible project-pitch figure set.

No output is a clinical score, diagnostic device, or treatment recommendation.

## Uncertainty treatment

Uncertainty is part of the result, not a final disclaimer. The project will separate:

- **sampling uncertainty:** small donor counts and within-donor cell sampling;
- **measurement uncertainty:** dropout, assay sensitivity, and tissue-processing effects;
- **mapping uncertainty:** ambiguous cell labels and out-of-reference disease states;
- **model uncertainty:** dependence on integration, graph, directionality, and prior choices;
- **transport uncertainty:** the distance between in vitro, blood, healthy tissue, and SCA tissue evidence;
- **missingness:** unavailable genotype, treatment, crisis-state, or phenotype variables.

Results should use intervals, donor-level effective sample sizes, calibration checks, and explicit sensitivity ranges. The system must be allowed to abstain. “No resolved estimate” is preferable to a precise number produced by an unidentifiable comparison.

## Ethics and fairness

SCA research has often been shaped by unequal access to care, uneven geographic sampling, and the misuse of race as a biological shortcut. This project will treat genotype, geography, environment, treatment access, and self-identified demographic variables as distinct concepts. Race will not stand in for genotype or ancestry, and ancestry will not be treated as destiny.

The dataset registry will expose representation by age, sex, genotype, geography, treatment, and clinical state when permitted. Performance and uncertainty will be reported by meaningful subgroups, without publishing unstable subgroup claims from inadequate samples. Controlled data will remain controlled; no attempt will be made to re-identify donors or infer sensitive attributes beyond consented use.

People living with SCA should be involved in deciding which outputs are useful and how disease states are described. Language should not depict a healthy reference as a universal human norm or a person with SCA as biologically reducible to disease. V1 is for research prioritization, not individual risk prediction, insurance, employment, or clinical triage.

## Initial MVP

The smallest credible demonstration will:

1. define and validate common contracts for source identity, donor/sample/cohort context, cell and anatomical semantics, provenance, access policy, and Q0–Q4 eligibility;
2. register one healthy adult heart reference, one healthy lung reference, GSE274340, GSE53441, and GSE84632 with immutable source identities and explicit evidence roles;
3. produce modality-appropriate Q3 views and declare which subsets, if any, reach Q4 for each proposed analysis;
4. expose one replayable cohort query and evidence bundle that joins permitted healthy reference, endothelial perturbation, and immune evidence without crossing access or consent boundaries;
5. focus the first fate-model demonstration on endothelial states and monocyte/macrophage programs, with mural and organ-resident cells providing signaling context;
6. estimate a small set of predeclared transition and communication probability shifts, validate at donor level, publish uncertainty and gap maps, and abstain where evidence is inadequate.

The MVP succeeds if a second analyst can reproduce the permitted evidence bundle and its traceable end-to-end inference, and if the exercise clearly identifies what matched SCA tissue data would be needed next. It need not discover a novel mechanism.

## Explicit non-goals

- Proving that SCA causes a modeled transition or signaling edge.
- Inferring literal lineage histories from cross-sectional transcriptomes.
- Building a whole-body SCA digital twin or simulating an individual patient.
- Diagnosing disease, predicting personal outcomes, or recommending treatment.
- Pooling all SCD genotypes, ages, therapies, and clinical states into one label.
- Generating new human data or performing wet-lab validation in v1.
- Exhaustively modeling every cardiac, pulmonary, blood, or immune population.
- Treating data colocation, format conversion, or a large integrated embedding as proof of interoperability.
- Exporting controlled records or derivatives beyond their consent and data-use boundaries.
- Selecting a production software stack, training a large foundation model, or launching a public portal before feasibility is established.

## Success criteria

The project advances beyond the MVP only if:

- every input passes a provenance and data-use audit, and every reported result traces back to immutable source files, semantic mappings, quality decisions, and an approved evidence bundle;
- repeated cohort queries produce equivalent permitted inputs or a complete, reviewable explanation of change;
- semantic and modality contracts reject incompatible or ambiguous inputs instead of silently coercing them;
- controlled-access analyses stay within their approved environments and exported derivatives pass the recorded policy rules;
- Q4 eligibility is analysis-specific, reproducible, and never inferred from study reputation or cell count alone;
- biological replication and validation operate at donor level;
- the main endothelial and immune findings retain their direction under leave-one-donor-out analysis and at least two reasonable integration or graph specifications;
- at least one prespecified program has compatible support in an independent SCA cohort, while contradictory evidence is reported plainly;
- disease-associated estimates are not presented when condition is inseparable from study, assay, tissue, or an unmeasured major covariate;
- subgroup coverage, mapping confidence, transport limitations, and model uncertainty appear beside the findings rather than in an appendix;
- domain and patient reviewers judge the claims understandable, non-stigmatizing, and no stronger than the evidence.

A scientifically useful null result—showing that available public data cannot identify a disease-conditioned transition—also meets the integrity goal if it produces a precise acquisition plan for the missing data.

## Proposed repository shape for later implementation

This is a proposal, not a request to create the directories yet.

```text
contracts/         Semantic, provenance, quality, policy, and interface schemas
configs/           Versioned cohort, ontology, and analysis settings
data/              Local or secure raw and derived data; never committed directly
docs/              Specifications, dataset/interface cards, and decisions
metadata/          Source registry, context records, and semantic crosswalks
notebooks/         Audits and bounded exploratory work
src/               Data-plane adapters and downstream model/reporting code
tests/             Contract, access-boundary, and analysis-invariant checks
results/           Permitted reproducible tables, figures, and result cards
```

The first implementation decision should be a review of the semantic and policy contracts against the candidate registry—not a storage engine or framework choice.
