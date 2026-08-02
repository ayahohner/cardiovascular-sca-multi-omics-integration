# Data

Source datasets and donor-level derived objects live here and are ignored by Git. The endothelial module freezes source metadata in `config/sources.json`. The HSPC fate module keeps source decisions in `research/hspc-fate/sources.json` and exact required-asset sizes and checksums in `research/hspc-fate/analysis.json`.

```bash
make fetch
make preflight
make hspc-fetch
make hspc-preflight
```

The admitted endothelial assets occupy about 6.5 GB. HSPC archives and rebuildable intermediates belong only under `data/raw/hspc-fate/` and `data/derived/hspc-fate/`. No raw or donor-level data belong in Git.
