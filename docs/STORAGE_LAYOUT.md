# Storage Layout

## Local Raw Archive

`PhaSePred_article&data/` contains local paper materials and extracted article
content. `data/raw/external/` contains local snapshots of public and licensed
data used for PhaSePred reconstruction. `tools/per-tool/` contains local external
feature-tool archives, source mirrors, and build outputs. These files are
intentionally ignored by git except for small README/manifest files.

Current local sources:

- PhaSePred paper PDFs, SI appendix, and supplementary tables S1-S8.
- PhaSepDB 3.0 API archive.
- PhaSePro full downloads.
- LLPSDB v2 zip packages and extracted tables.
- DeepPhase GitHub supplement and extracted tables.
- PhosphoSitePlus phosphorylation dataset.
- External feature-tool material under `tools/per-tool/`, including PScore, CD-HIT,
  PLAAC, ESpritz, IUPred3, catGRANULE v2, and the NCBI SEG source mirror.

See `docs/DATA_SOURCES.md` for details.

## Checksums

The local checksum file is:

```bash
data/raw/external/_manifests/sha256sums.txt
```

It may include ignored large files, so it remains local by default. Regenerate it
after refreshing data:

```bash
find data/raw/external -type f | sort | xargs -r sha256sum > data/raw/external/_manifests/sha256sums.txt
```

## Future Derived Data

Use:

- `data/interim/` for deterministic intermediate tables.
- `data/processed/` for final train/test tables.
- `models/` for model artifacts if local storage is acceptable.
- `runs/` or `outputs/` for experiment logs and predictions.

These paths are ignored until a storage policy says otherwise.

The 2026-05-11 ESpritz-vs-IUPred3 IDR comparison uses:

- `data/interim/uniprot_cache.jsonl` for UniProt sequence records.
- `data/interim/iupred3_idr_cache.jsonl` for IUPred3-derived disorder
  fractions.
- `runs/iupred_idr_comparison/` for models, predictions, AUC tables, missing
  sequence logs, and figures.

These are reproducible generated artifacts and should remain out of normal git
history. The reproduction command and AUC table are documented in
`docs/EXPERIMENTS.md`.

The 2026-05-11 external catGRANULE2 and DeepCoil scoring work uses:

- `data/interim/catgranule2_human_features.csv` for the bundled catGRANULE2
  human-proteome scores mapped to recomputed sequence accessions.
- `data/interim/catgranule2_human_missing.csv` for accessions absent from that
  bundled human-proteome table.
- `data/interim/deepcoil_features.csv` for any DeepCoil aggregate table
  produced by a batch run.
- `runs/deepcoil_scores/` for resumable DeepCoil batch FASTA files, raw output,
  and per-batch logs.
- `.external_envs/catgranule2` and `.external_envs/deepcoil` for isolated conda
  prefixes. These are local runtime artifacts and should not be committed.
