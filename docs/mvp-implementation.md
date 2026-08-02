# MVP implementation and traceability

**Release:** `mvp-0.1.0` · **Pipeline:** `0.1.0` · **Scope:** focused open-data phase

## Executed slice

The MVP executes one complete scientific increment:

1. admit exact processed assets from the Heart Cell Atlas vascular object, healthy HLCA core, GEO GSE274340, GEO GSE53441, and the GPL570 probe annotation;
2. aggregate raw counts into donor–state pseudobulks, normalize to log-CPM, and retain organ, state, donor, sample, study, and source context;
3. define organ-specific mutual-nearest-neighbor adjacencies on 1,500 state-variable genes after excluding mitochondrial and ribosomal prefixes;
4. require at least three donors, cosine distance at most 0.8, donor leave-out stability of at least 0.6, study leave-out stability of at least 0.5 where applicable, and parameter stability of at least 0.5;
5. derive a descriptive HUVEC perturbation signature from GSE274340 and a donor-sample PBMC association signature from GSE53441;
6. map the two evidence roles separately to healthy state contrasts using absolute cosine alignment; and
7. rank only eligible adjacencies with visible component scores, uncertainty, counterevidence, claim boundaries, and next-data requirements.

The state-A-to-state-B contrast is an arbitrary record orientation. Its signed alignment is not temporal direction. All released edges are explicitly undirected.

## Requirement traceability

| EARS requirements | MVP status | Executable evidence |
| --- | --- | --- |
| REQ-001–006 | Implemented | `config/analysis.json`, state reference, scope and claim checks |
| REQ-007–012 | Implemented | `config/sources.json`, admission preflight, source and provenance manifests, checksum quarantine tests |
| REQ-013–016 | Implemented with one abstention | 14 supported states; no defensible endocardial source label, recorded as an abstention; 28-cell proliferating lung lymphatic state excluded |
| REQ-017–022 | Implemented | adjacency gates, undirected evidence fields, separate HUVEC/PBMC roles, overclaim failure tests |
| REQ-023–025 | Implemented for the local MVP | sparse backed reads, 2,048-cell batches, checksummed caches, measured memory/runtime/disk preflight |
| REQ-026–027 | Implemented | donor and study leave-outs, two feature-set sizes, three neighborhood sizes, explicit exclusions |
| REQ-028–034 | Implemented | component rankings, uncertainty, decision log, claim-boundary checks, release and follow-up manifests |

REQ-003 and REQ-019 are optional evidence paths. No adult-organ public velocity, time-series, or direct perturbation object was admitted in the MVP, so the correct output is an undirected edge with absent directional evidence.

## Actual controls and result

- Five source assets pass exact size and SHA-256 verification.
- Fourteen states pass donor–state eligibility: seven heart and seven lung.
- One expected heart state is an abstention because the source labels do not support an endocardial mapping.
- Ten mutual-neighbor edges were evaluated. Nine passed; the cardiac atrial–lymphatic edge failed the distance gate (`0.912 > 0.8`).
- Both SCA evidence roles have hundreds of overlapping genes for each tested contrast, but remain indirect transport evidence.
- Fifteen tests pass, including two complete real-data releases with identical scientific artifact hashes. Coverage is 87.5% against an 85% gate.
- The release verifier checks every compact artifact by size and SHA-256 and confirms that raw and donor-level data are not tracked.

## Remaining evidence gaps

The MVP does not supply direction, matched SCA cardiopulmonary tissue, longitudinal observation, lineage evidence, or direct disease-effect estimation. GSE274340 has one sample per condition and uses a HUVEC microfluidic model. GSE53441 is bulk PBMC data and contains a documented discrepancy between reported healthy-subject count and control-sample count. These limitations lower transport confidence and remain visible in every ranked result.

The next useful increment is not a broader platform. It is an independent, donor-resolved validation source for one leading organ-specific edge, followed by directionality only if an eligible time-resolved or direct perturbation design exists.
