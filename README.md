# Cardiovascular SCA Multi-Omics Integration

A working open-data synthesis for prioritizing endothelial state-transition hypotheses relevant to cardiovascular sickle cell anemia (SCA) research.

## Current phase

The MVP now:

1. build a healthy heart-and-lung endothelial state reference from open processed public data;
2. define reproducible, biologically adjacent candidate state transitions;
3. add public directionality or perturbation evidence only where the data support it;
4. map public SCA-relevant endothelial and immune signals onto those transitions; and
5. rank uncertainty-qualified hypotheses and specify the data needed to test them.

The first release assembled 14 supported healthy endothelial states, evaluated 10 candidate adjacencies, rejected one distant forced-neighbor artifact, and retained 9 uncertainty-qualified hypotheses. Every edge remains undirected. The leading results are prioritization targets, not observed fate or evidence that SCA changes a transition.

This phase does **not** estimate fate in patient SCA tissue or claim that SCA causes transition shifts. Matched-SCA-tissue testing, longitudinal or lineage proof, generalized platform work, controlled-access data, ligand–receptor networks, and clinical prediction remain deferred.

See the [MVP result](results/mvp/report.md), [implementation and traceability notes](docs/mvp-implementation.md), [scientific specification](docs/disease-conditioned-cellular-fate-model.md), and [EARS requirements](docs/current-phase-ears-requirements.md).

## Reproduce

Python 3.12 and about 7 GB of source-data storage are required.

```bash
make setup
make fetch
make preflight
make mvp
make verify
make test
```

`make fetch` retrieves five pinned public assets and verifies their byte counts and SHA-256 checksums. Raw and donor-level derived data remain outside Git. The compact release in `results/mvp/` contains source and provenance manifests, state centroids and programs, candidate adjacencies, SCA evidence mappings, component-level rankings, control decisions, validation evidence, and the follow-up plan.

Use `make test-fast` for tests that do not require the public assets. The full suite runs the real pipeline twice and requires scientific artifact hashes to match; it currently passes 15 tests with 87.5% coverage.

## Project history

The workspace began with spatially constrained cardiomyocyte–fibroblast signaling after myocardial infarction, with emphasis on how local tissue context shapes intercellular communication. The current project direction is the focused SCA-relevant endothelial analysis described above.

## Workspace layout

```text
data/       Local inputs and intermediate data (not committed)
config/     Frozen source registry and analysis thresholds
docs/       Scientific, system, and implementation specifications
src/        Reusable Python pipeline
tests/      Unit, control-failure, and real-data integration tests
results/    Compact, versioned MVP release
```
