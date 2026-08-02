# Data

Source datasets and donor-level derived objects live here and are ignored by Git. The implemented endothelial module freezes its source metadata in `config/sources.json`. The planned HSPC fate module keeps pointers and admission decisions in `research/hspc-fate/sources.json`; retrieval will freeze exact file checksums before analysis.

```bash
make fetch
make preflight
```

The admitted endothelial assets occupy about 6.5 GB. Their existing paths under `data/raw/` and `data/derived/` remain unchanged. Future HSPC inputs and rebuildable outputs belong only under `data/raw/hspc-fate/` and `data/derived/hspc-fate/`. No raw or donor-level data belong in Git.
