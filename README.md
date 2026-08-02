# Cardiovascular SCA Multi-Omics Integration

A focused open-data synthesis project for prioritizing endothelial state-transition hypotheses relevant to cardiovascular sickle cell anemia (SCA) research.

## Current phase

The current phase will:

1. build a healthy heart-and-lung endothelial state reference from open processed public data;
2. define reproducible, biologically adjacent candidate state transitions;
3. add public directionality or perturbation evidence only where the data support it;
4. map public SCA-relevant endothelial and immune signals onto those transitions; and
5. rank uncertainty-qualified hypotheses and specify the data needed to test them.

This phase does **not** estimate fate in patient SCA tissue or claim that SCA causes transition shifts. Matched-SCA-tissue testing, longitudinal or lineage proof, generalized platform work, controlled-access data, ligand–receptor networks, and clinical prediction are deferred.

See the [current-phase project specification](docs/disease-conditioned-cellular-fate-model.md) for the evidence rules, failure controls, outputs, and success criteria. Deferred areas are not part of implementation.

## Project history

The workspace began with spatially constrained cardiomyocyte–fibroblast signaling after myocardial infarction, with emphasis on how local tissue context shapes intercellular communication. The current project direction is defined in the specification linked above; no analysis stack has been selected yet.

## Workspace layout

```text
data/       Local inputs and intermediate data (not committed)
notebooks/  Exploratory analyses
src/        Reusable analysis code
results/    Generated figures, tables, and other outputs (not committed)
```

No language, package manager, or analysis framework has been selected yet.
