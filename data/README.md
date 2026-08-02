# Data

Source datasets and donor-level derived objects live here and are ignored by Git. Their immutable accessions, URLs, sizes, SHA-256 checksums, terms, citations, covariates, and limitations are frozen in `config/sources.json`.

```bash
make fetch
make preflight
```

The admitted assets occupy about 6.5 GB. `data/raw/` contains the five verified source files. `data/derived/` contains rebuildable pseudobulk reference caches. No raw or donor-level data belong in Git; the compact state-level release is under `results/mvp/`.
