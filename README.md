# Cardiovascular SCA Multi-Omics Integration

Open-data research modules for cardiovascular sickle cell anemia (SCA). Each module owns its question, sources, methods, and results; modules do not depend on one another.

## Research modules

- **Endothelial transitions — implemented.** The released MVP prioritizes healthy heart-and-lung endothelial state adjacencies using indirect SCA-relevant evidence. See the [result](results/mvp/report.md), [implementation notes](docs/mvp-implementation.md), [scientific specification](docs/disease-conditioned-cellular-fate-model.md), and [EARS requirements](docs/current-phase-ears-requirements.md).
- **HSPC lineage priming — implemented v2.** The independent study tests donor-level SCD associations against a donor-balanced healthy reference using Palantir, CellRank-CytoTRACE, OOD rejection, held-out markers, and full model refits. Its valid result remains `insufficient_evidence`. See the [full study report](research/hspc-fate/study-report.md), [MVP-to-v2 migration](research/hspc-fate/migration-v2.md), [compact result](results/hspc-fate/report.md), and [module contract](research/hspc-fate/README.md).

The [research index](research/README.md) defines the module boundary. The existing endothelial implementation remains the default executable project.

## Endothelial MVP

The MVP now:

1. build a healthy heart-and-lung endothelial state reference from open processed public data;
2. define reproducible, biologically adjacent candidate state transitions;
3. add public directionality or perturbation evidence only where the data support it;
4. map public SCA-relevant endothelial and immune signals onto those transitions; and
5. rank uncertainty-qualified hypotheses and specify the data needed to test them.

The first release assembled 14 supported healthy endothelial states, evaluated 10 candidate adjacencies, rejected one distant forced-neighbor artifact, and retained 9 uncertainty-qualified hypotheses. Every edge remains undirected. The leading results are prioritization targets, not observed fate or evidence that SCA changes a transition.

This phase does **not** estimate fate in patient SCA tissue or claim that SCA causes transition shifts. Matched-SCA-tissue testing, longitudinal or lineage proof, generalized platform work, controlled-access data, ligand–receptor networks, and clinical prediction remain deferred.

## Reproduce

Python 3.12 and about 9 GB of source-data storage are required to reproduce both modules.

```bash
make setup
make fetch
make preflight
make mvp
make verify
make hspc-fetch
make hspc-preflight
make hspc-mvp
make hspc-verify
make test
```

`make fetch` retrieves the endothelial assets; `make hspc-fetch` retrieves seven donor-separated or donor-resolved HSPC archives. Both workflows verify frozen byte counts and SHA-256 checksums. Raw and donor-level derived data remain outside Git. Compact, auditable releases live in `results/mvp/` and `results/hspc-fate/`.

Use `make test-fast` for tests that do not require the public assets. The full suite reruns the real pipelines, checks source formats, QC and donor balancing, OOD mapping, donor-level inference, model robustness records, and deliberate release corruption.

## Project history

The workspace began with spatially constrained cardiomyocyte–fibroblast signaling after myocardial infarction. The endothelial synthesis and HSPC fate study remain intentionally separate scientific increments.

## Workspace layout

```text
data/       Local inputs and intermediate data (not committed)
config/     Frozen source registry and analysis thresholds
docs/       Scientific, system, and implementation specifications
research/   Independent research-module contracts and source pointers
src/        Reusable Python pipeline
tests/      Unit, control-failure, and real-data integration tests
results/    Compact, versioned MVP release
```
