# Open-Data HSPC Fate Study: v1 Report

- **Project:** Cardiovascular SCA Multi-Omics Integration
- **Module:** HSPC fate
- **Release:** v1
- **Analysis date:** 2 August 2026
- **Released conclusion:** `insufficient_evidence`

## Executive summary

This study asked whether open human single-cell data support a reproducible association between sickle cell disease and the inferred tendency of hematopoietic stem and multipotent progenitor cells (HSCs/MPPs) to allocate toward three broad blood-cell lineages: myeloid, erythroid–megakaryocytic, and lymphoid.

The workflow was completed successfully. Seven public archives were downloaded and verified, 103,095 quality-controlled cells were analyzed across two independent cohorts, healthy-only cell-state models passed their topology controls, SCD cells were mapped without changing those models, and all donor-level estimates and uncertainty checks were released. The final scientific conclusion is nevertheless `insufficient_evidence`.

That conclusion has two concrete reasons:

1. None of the three discovery-cohort lineage contrasts passed every predeclared stability control.
2. In the independent validation cohort, donor SCD5 retained only 4 HSC/MPP-like cells, below the frozen reliability minimum of 20. Its estimate is available for audit but is too fragile to support replication.

This is a useful result, not a failed analysis. It establishes a reproducible open-data route to disease-conditioned HSPC allocation estimates, shows exactly where the present evidence stops, and identifies the smallest data acquisition needed next. It does **not** show observed cell fate, prove that sickle cell disease causes a transition, or support clinical prediction.

## 1. Why this study was undertaken

The broader repository studies open-data questions relevant to sickle cell anemia (SCA). Public cardiovascular tissue data were not sufficient for direct SCA fate estimation, so this module pursued a narrower question for which donor-resolved human SCD data actually exist.

HSPCs are the cells that maintain blood production. An HSC/MPP can be described as being closer to several possible lineage programs, but a single cross-sectional RNA profile does not reveal which lineage that cell later entered. The term **fate probability** in this report therefore means a model-based allocation across reference endpoints. It is not an observed future outcome.

This HSPC module is scientifically independent from the repository's endothelial-transition module. The two modules do not share data, state definitions, estimands, model fitting, or release gates. HSPCs are not cardiovascular tissue, and their results must not be presented as direct evidence about cardiovascular cell fate.

## 2. Research question and estimand

The precise question was:

> Among human HSC/MPP-like cells, is sickle cell disease associated with a reproducible difference in inferred allocation toward myeloid, erythroid–megakaryocytic, or lymphoid progenitor fates?

For each donor and lineage, the workflow calculated the mean inferred terminal-fate probability among that donor's HSC/MPP-like cells. The condition effect was then:

`mean of SCD donor estimates − mean of healthy donor estimates`

Positive values mean higher inferred allocation in the SCD group; negative values mean lower inferred allocation. Donors—not cells—were the biological replicates. This prevents thousands of cells from being incorrectly treated as thousands of independent people.

The source datasets use the broader label SCD, or sickle cell disease. The repository name uses SCA, or sickle cell anemia. Condition labels in the analysis retain the sources' `SCD` terminology.

## 3. Scope

### Included

- Open human single-cell transcriptomic data with donor and condition labels.
- A discovery comparison with 3 healthy and 3 SCD donors.
- An independent validation comparison with 4 healthy and 5 SCD donors.
- Healthy-only reference models for HSC/MPP and three terminal-lineage families.
- Two independently structured allocation estimators.
- Donor-level contrasts, exact label permutations, leave-one-donor-out analysis, donor bootstrap intervals, and explicit failure gates.

### Excluded

- Observed or longitudinal cell fate.
- Lineage tracing and RNA velocity without compatible molecular layers.
- Causal claims about SCD.
- Cardiovascular tissue inference.
- Treatment response, gene-therapy efficacy, individual prediction, or clinical use.
- Controlled-access datasets.
- Pooling discovery and validation to manufacture a larger apparent sample.
- A generalized multi-omics platform.

## 4. Public data

All required data were accessible without an institutional agreement or research budget when admitted on 2 August 2026.

| Role | Accession | Material | Included comparison | Access and use |
|---|---|---|---|---|
| Discovery | [S-BSST1258](https://www.ebi.ac.uk/biostudies/studies/S-BSST1258) | Donor-separated, unfiltered 10x scRNA-seq matrices from plerixafor-mobilized CD34-positive HSPCs | 3 healthy and 3 SCD donors | Open; BioStudies reports CC0 |
| Validation | [GSE133181](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133181) | Donor-resolved 10x scRNA-seq matrices from bone-marrow CD34-positive cells | 4 healthy and 5 SCD donors | Open under NCBI GEO repository terms |

Three thalassemia samples and one lineage-depleted healthy sample in GSE133181 were excluded before analysis because they do not answer the healthy-versus-SCD estimand.

The discovery dataset accompanies *Severe inflammation and lineage skewing is associated with poor engraftment of engineered hematopoietic stem cells in patients with sickle cell disease* ([Nature Communications, 2025](https://doi.org/10.1038/s41467-025-58321-4)). The validation dataset accompanies Hua et al., *Single-cell RNA-seq analysis of human bone marrow reveals molecular features of hematopoietic stem cells in sickle cell disease* ([Blood, 2020](https://doi.org/10.1182/blood.2019002301)).

Two additional open sources were retained as future pointers but were not used to produce v1 fate estimates:

- [GSE277682](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277682), human bulk HSPC RNA-seq, can qualify program direction but cannot estimate cell-state allocation.
- [GSE277681](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277681), mouse Townes-model HSPC CITE-seq, can provide model support but is not human replication.

Longitudinal gene-therapy clonal-tracking data were excluded because the raw sequencing is controlled, the public derived repository has no explicit reuse license, and the study lacks the healthy comparator required by this estimand.

### Provenance and integrity

The six discovery ZIP files and one validation TAR file total approximately 1.9 GiB. Every required archive has a frozen resolved URL, expected byte count, SHA-256 checksum, retrieval time, terms, donor label, assay, tissue, preprocessing state, and limitations. All seven passed byte-for-byte admission before analysis.

Raw archives and donor-level derived matrices remain outside Git. The compact release contains the information needed to identify and verify them. See the [source manifest](../../results/hspc-fate/source_manifest.json), [donor manifest](../../results/hspc-fate/donor_manifest.tsv), and [frozen analysis configuration](analysis.json).

## 5. Analysis method

### 5.1 Donor-resolved preprocessing

Each donor was processed separately. The discovery archives contain unfiltered feature–barcode matrices, so the workflow applied:

- at least 3,000 total UMIs;
- at least 500 detected genes;
- no more than 15% mitochondrial counts;
- detectable `CD34` expression; and
- Scrublet doublet removal.

A barcode-rank knee was recorded for each donor as a diagnostic. The validation matrices were already source-filtered, so they used a 1,000-UMI floor, at least 500 genes, no more than 15% mitochondrial counts, and Scrublet. CD34 expression was not reapplied as a hard validation filter.

This is deliberately not described as an exact reproduction of either source paper. The discovery workflow uses the paper's published 3,000-UMI floor but does not rerun EmptyDrops at its reported false-discovery threshold. Scrublet replaces the paper's DoubletFinder implementation. These deviations are visible in the release rather than hidden.

| Cohort | Donor | Condition | QC-retained cells | HSC/MPP-like cells | Root-cell gate |
|---|---:|---|---:|---:|---|
| Discovery | HD175 | Healthy | 8,244 | 2,824 | Pass |
| Discovery | HD177 | Healthy | 9,755 | 2,379 | Pass |
| Discovery | HD192 | Healthy | 10,845 | 3,526 | Pass |
| Discovery | SCD2 | SCD | 8,038 | 3,561 | Pass |
| Discovery | SCD3 | SCD | 11,047 | 4,816 | Pass |
| Discovery | SCD4 | SCD | 11,042 | 2,970 | Pass |
| Validation | BM1 | Healthy | 10,757 | 4,067 | Pass |
| Validation | BM2 | Healthy | 11,166 | 3,310 | Pass |
| Validation | BM3 | Healthy | 5,005 | 1,743 | Pass |
| Validation | BM4 | Healthy | 5,094 | 644 | Pass |
| Validation | SCD1 | SCD | 2,311 | 103 | Pass |
| Validation | SCD2 | SCD | 2,282 | 356 | Pass |
| Validation | SCD3 | SCD | 3,142 | 672 | Pass |
| Validation | SCD4 | SCD | 4,207 | 279 | Pass |
| Validation | SCD5 | SCD | 160 | 4 | **Fail** |

The discovery cohort retained 58,971 cells, 7.8% above the source publication's reported 54,689. Validation retained 44,124 cells, 8.0% above the publication's reported 40,862. Both differences are inside the frozen 25% source-reproduction tolerance. Differences are expected because the lean workflow is reproducible but not method-identical to the publications.

One Scrublet fit for validation donor BM3 emitted a numerical runtime warning while still returning finite scores and calls. The workflow captured the warning and recorded the threshold mode rather than suppressing it silently. The full donor preprocessing and release build completed successfully.

### 5.2 Healthy-only state references

Discovery and validation were modeled independently. In each cohort, only healthy cells were allowed to define the reference geometry and endpoint semantics. SCD cells could be mapped into the reference but could not move the reference.

The model used:

- log-normalized expression;
- 1,000 highly variable non-mitochondrial, non-ribosomal genes;
- 25 principal components;
- 30 healthy reference states; and
- predeclared marker programs for an HSC/MPP-like root and myeloid, erythroid–megakaryocytic, and lymphoid endpoints.

Examples of root markers include `HLF`, `MECOM`, `MLLT3`, `AVP`, and `GATA2`. Examples of terminal markers include `MPO` and `ELANE` for myeloid, `GATA1` and `KLF1` for erythroid–megakaryocytic, and `DNTT`, `IL7R`, and `EBF1` for lymphoid. The complete frozen marker sets are in [analysis.json](analysis.json).

Both healthy-only topologies passed their required controls:

| Control | Discovery | Validation |
|---|---:|---:|
| Healthy donors supporting the root | 3 | 4 |
| Healthy donors supporting myeloid endpoint | 3 | 3 |
| Healthy donors supporting erythroid–megakaryocytic endpoint | 3 | 2 |
| Healthy donors supporting lymphoid endpoint | 3 | 2 |
| Root reaches all three terminal families | Pass | Pass |
| Terminal markers increase relative to root | Pass | Pass |
| Neighbor-sensitivity correlation | 0.972 | 0.984 |
| Required neighbor-sensitivity correlation | 0.800 | 0.800 |

Hashes of the transformation, topology, root definition, terminal definitions, and model parameters were identical before and after SCD mapping in both cohorts. This verifies that disease-labeled cells did not refit the healthy model.

### 5.3 Two allocation estimators

The same cells were evaluated in two different ways:

1. **Absorbing Markov estimator.** A directed healthy state graph moves from earlier to later states and calculates each state's probability of eventually reaching the three terminal families.
2. **Terminal projection estimator.** Cells are projected toward pseudotime-defined healthy terminal centroids without using the Markov graph.

Both estimators were required to return three finite, nonnegative probabilities that sum to one. Agreement between their effect directions was used as a guard against a result being created by one modeling choice.

### 5.4 Donor-level uncertainty and decision rules

Cell probabilities were averaged within each donor before comparing conditions. The workflow then recorded:

- all 20 possible label assignments in the 3-versus-3 discovery cohort;
- all 126 possible label assignments in the 4-versus-5 validation cohort;
- the effect after leaving out each donor in turn;
- agreement between the two estimators; and
- 2,000 donor-bootstrap replicates for a descriptive 95% interval.

These controls assess stability; they do not create a large-sample significance claim. The exact-permutation tail fraction is especially coarse with so few donors. A discovery lineage was considered stable only when estimator direction agreed, at least 5 of 6 leave-one-donor-out effects retained the direction, and both frozen permutation screens passed.

## 6. Results

### 6.1 Discovery cohort

Effects below are SCD minus healthy and are shown in probability points. For example, `−0.025` is approximately −2.5 percentage points of model allocation.

| Lineage | Markov effect | Projection effect | Estimator direction | Permutation tail, M/P | Leave-one-out direction, M/P | Result |
|---|---:|---:|---|---:|---:|---|
| Myeloid | −0.0250 | −0.0047 | Agree: lower | 0.70 / 0.30 | 5/6 / 6/6 | Unstable |
| Erythroid–megakaryocytic | +0.0230 | +0.0089 | Agree: higher | 0.70 / 0.50 | 5/6 / 6/6 | Unstable |
| Lymphoid | +0.0021 | −0.0043 | Disagree | 0.90 / 0.50 | 4/6 / 5/6 | Unstable |

All three Markov bootstrap intervals crossed zero:

- myeloid: −0.0929 to +0.0641;
- erythroid–megakaryocytic: −0.0763 to +0.0892; and
- lymphoid: −0.0108 to +0.0104.

The myeloid and erythroid–megakaryocytic estimators agreed on direction, but label permutations produced effects at least as strong too often for those contrasts to pass the frozen null screen. The lymphoid effect was small, changed direction between estimators, and depended more strongly on individual donors. Therefore, no discovery lineage was stable.

### 6.2 Validation cohort

The validation topology recovered the same root and terminal families, but donor SCD5 failed the minimum root-cell gate. The following numbers are retained for transparency and method checking only:

| Lineage | Markov effect | Projection effect | Estimator direction | Internal stability calculation |
|---|---:|---:|---|---|
| Myeloid | −0.0337 | +0.0066 | Disagree | Unstable |
| Erythroid–megakaryocytic | +0.0581 | +0.0116 | Agree | Passes internal calculation; exploratory only |
| Lymphoid | −0.0244 | −0.0181 | Agree | Passes internal calculation; exploratory only |

The erythroid–megakaryocytic and lymphoid validation patterns cannot be called replicated results. They were not preceded by stable discovery contrasts, and the validation donor-reliability gate failed. Removing SCD5 after seeing the result, pooling the cohorts, or lowering the 20-cell floor would change the frozen study contract and risk a biased conclusion.

### 6.3 Final classification

The workflow supports four possible conclusions: replicated difference, discovery-specific difference, no stable difference, or insufficient evidence. A failed required source, topology, annotation, or inference gate overrides a positive or null contrast classification.

Nine of ten release gates passed. The failed gate was:

`validation_root_cell_reliability: SCD5 retained 4 HSC/MPP-like cells; minimum = 20`

The deterministic v1 conclusion is therefore:

> **`insufficient_evidence`** — the available open data do not support a reliable replicated SCD-associated difference in inferred HSPC lineage allocation under the frozen v1 model and gates.

This is not equivalent to proving that no biological difference exists. It means that this design and these public cohorts cannot distinguish a reproducible difference with adequate reliability.

## 7. What the result does and does not mean

### Supported interpretation

The release describes cross-sectional donor-level associations in model-inferred lineage allocation among transcriptomically HSC/MPP-like cells. It also demonstrates that a complete, reproducible, zero-budget analysis is technically feasible with open data.

### Unsupported interpretations

The analysis does not establish:

- that any individual cell entered a particular lineage;
- that SCD caused a change in lineage transition probability;
- that a lineage difference is replicated;
- that mobilization, treatment, age, tissue source, or other cohort contexts are irrelevant;
- that an estimate applies to an individual patient; or
- that the model has clinical or treatment-selection value.

## 8. Principal limitations

1. **Small donor counts.** Six discovery and nine validation donors produce coarse exact-permutation resolution and wide donor-bootstrap intervals.
2. **A failed validation donor gate.** Four root-like cells cannot provide a reliable donor mean for SCD5.
3. **Different collection contexts.** Discovery uses plerixafor-mobilized HSPCs collected in a gene-therapy context; validation uses pediatric bone-marrow CD34-positive cells from a distinct study.
4. **Cross-sectional transcriptomics.** The models infer proximity to endpoint programs; they do not observe a transition through time.
5. **Method deviations.** The lean QC implementation is transparent and reproducible but not identical to either source publication.
6. **Incomplete cohort metadata.** Ancestry, geography, environment, access to care, and self-identified race are unavailable at the required donor level and remain `not_reported`.
7. **Broad endpoint families.** Three lineage families provide a tractable v1 question but do not represent every HSPC state or differentiation path.

## 9. Ethics and fairness

SCD genotype is not used as a proxy for ancestry, race, geography, environment, or access to care. The analysis does not infer missing demographic attributes, and it does not release donor-level sensitive metadata beyond the public sources. Differences between cohorts are treated as transport limitations, not as nuisance variation to be pooled away.

Because collection and treatment contexts may correlate with condition, even a stable association would not establish a disease mechanism. Any future causal claim would require a design with stronger temporal, experimental, or matched-tissue evidence.

## 10. Reproducibility and software validation

The HSPC module is namespaced under `research/hspc-fate/`, `src/hspc_fate_mvp/`, `data/*/hspc-fate/`, and `results/hspc-fate/`. It does not depend on the endothelial module's scientific artifacts.

Reproduction requires Python 3.12, standard open-source Python packages, internet access for the public archives, and approximately 1.9 GiB for HSPC source downloads:

```bash
make setup
make hspc-fetch
make hspc-preflight
make hspc-mvp
make hspc-verify
make test
```

Validation at release included 34 tests with 85.48% combined branch coverage. The tests exercise archive parsing, source admission, QC exclusion accounting, doublet fallback behavior, graph absorption, donor-level aggregation, exact permutations, classification, command failure exits, release integrity, deliberate probability corruption, and two identical real-data runs. The repeated runs produced identical scientific artifact hashes.

The [release manifest](../../results/hspc-fate/release_manifest.json) binds every compact artifact to its byte count and SHA-256 checksum. The independent release verifier passed with no problems.

## 11. Recommended next step

The smallest scientifically useful next increment is not a larger platform or a more elaborate model. It is one additional open, donor-resolved human SCD HSPC cohort that can satisfy the existing reliability contract.

The sequence should be:

1. Identify an independent open human HSPC dataset with healthy and SCD donors.
2. Verify terms, donor labels, context, archive integrity, and reuse boundaries before analysis.
3. Apply the frozen lineage vocabulary and estimate how many HSC/MPP-like cells remain per donor.
4. Continue only if every comparator donor has at least 20 root-like cells.
5. Re-run the same donor-level estimand and both allocation estimators without changing endpoints or pooling cohorts.
6. Report concordance, disagreement, or another stopped result under the same claim boundary.

Human bulk RNA-seq or mouse CITE-seq may be added later to qualify lineage, stemness, or senescence programs. Neither can replace independent human single-cell replication or rescue the failed validation gate.

## 12. Artifact guide

| Artifact | Purpose |
|---|---|
| [Frozen analysis configuration](analysis.json) | Required archives, QC thresholds, marker programs, model parameters, and claim limits |
| [Formal v1 specification](v1-spec.md) | EARS requirements and acceptance criteria |
| [Operation sequence](operations.md) | Discrete workflow, pass conditions, and failure responses |
| [Source manifest](../../results/hspc-fate/source_manifest.json) | Resolved public sources, integrity, terms, and limitations |
| [Donor manifest](../../results/hspc-fate/donor_manifest.tsv) | Donor, condition, tissue, inclusion, and missing-context fields |
| [QC ledger](../../results/hspc-fate/qc_ledger.tsv) | Donor-level filtering and diagnostic details |
| [State definitions](../../results/hspc-fate/state_definitions.json) | Root, terminal families, states, and directed graphs |
| [Topology controls](../../results/hspc-fate/topology_controls.json) | Connectivity, donor support, sensitivity, and frozen model hashes |
| [Donor fate probabilities](../../results/hspc-fate/donor_fate_probabilities.tsv) | One estimate per donor, estimator, and lineage |
| [Stability results](../../results/hspc-fate/stability_results.tsv) | Effects, every exact permutation, leave-one-out effects, intervals, and decisions |
| [Replication table](../../results/hspc-fate/replication_table.tsv) | Separate discovery and validation comparison |
| [Gate ledger](../../results/hspc-fate/gate_ledger.tsv) | Passed and failed scientific controls |
| [Classification](../../results/hspc-fate/result_classification.json) | Machine-readable final conclusion |
| [Claim boundary](../../results/hspc-fate/claim_boundary.md) | Permitted and excluded interpretations |
| [Compact result summary](../../results/hspc-fate/report.md) | Short release-facing result |

## Conclusion

The project achieved its lean engineering goal: a reproducible, provenance-bound workflow that converts open human SCD HSPC data into donor-level, uncertainty-qualified allocation estimates. It also achieved the more important scientific goal of refusing to overstate those estimates.

The present evidence does not support a replicated SCD-associated HSPC fate-allocation difference. The discovery effects are unstable, and the independent validation comparison fails a donor-level reliability gate. The correct next move is to acquire one suitable open human cohort and rerun the frozen question—not to reinterpret exploratory validation patterns as proof.
