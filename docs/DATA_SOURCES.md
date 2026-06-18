# Data Sources

Generated source manifest:

- `data/raw/external/_manifests/RAW_DATA_MANIFEST.md`
- `data/raw/external/_manifests/sha256sums.txt`

`sha256sums.txt` is local-only by default because it references ignored raw
archive files. Regenerate it after refreshing local data.

## 2022 paper archive (this lab's prior work)

Path: `PhaSePred_article&data/`

Contains the published PDF, SI appendix PDF, and supplementary datasets S1-S8
from the Tingting Li Lab's 2022 PhaSePred paper (Chen et al., *PNAS*
119(24):e2115369119). These are enough to train a baseline from the already-
computed feature tables in S2/S3, but not enough to recover the original
trained model weights.

The directory is ignored by git except for its README. Keep the local files in
place when working on the 2022-protocol baseline (Product A).

## PhaSepDB 3.0

Path: `data/raw/external/phasepdb3/`

Archived from `https://db.phasep.pro/api`:

- `proteins.jsonl`: 3528 records from `/proteins`.
- `uniprot.jsonl`: 1873 records from `/uniprot`.
- `mlos.json`: MLO records from `/mlos`.
- `stats.json`: API stats snapshot.
- `manifest.json`: endpoint counts from the archive run.

Known count difference: the API stats snapshot reports 3742 proteins and 160
MLOs, while the archived endpoints returned 3528 `/proteins` records and 124
`/mlos` records. Keep this discrepancy visible in derived-data documentation.

## PhaSePro

Path: `data/raw/external/phasepro/`

Archived from official download endpoints:

- `download_full.json`
- `download_full.tsv`
- `download_full.xml`

## LLPSDB v2

Path: `data/raw/external/llpsdb2/`

Archived official zip packages:

- `Phase_separation_ambiguous.zip`
- `Phase_separation_unambiguous.zip`
- `No_phase_separation_ambiguous.zip`
- `No_phase_separation_unambiguous.zip`
- `Phase_diagram_ambiguous.zip`
- `Phase_diagram_unambiguous.zip`

Extracted files are under `data/raw/external/llpsdb2/extracted/`. Each package
contains `protein.xls` and `LLPS.xls`.

## DeepPhase

Path: `data/raw/external/deepphase/`

Archived from `https://github.com/cheneyyu/DeepPhase` at commit
`f19659134688542392c56dd45f0da915c67111a8`.

Files include:

- `DeepPhase_supp.zip`
- `README.md`
- `LICENSE`
- `dropcount.cpproj`
- GitHub metadata JSON files

The supplement zip is extracted under `data/raw/external/deepphase/extracted/`
and contains `tableS1.xlsx` through `tableS8.xlsx` plus supplemental figures.
`tableS3.xlsx`, sheet `tableS3`, contains `Swiss-Prot ID` and
`DeepPhase_score`; the `features-recomputed` workflow maps those columns to
`UniprotEntry` and `DeepPhase`.

## PhosphoSitePlus

Path: `data/raw/external/phosphositeplus/`

Downloaded with user permission:

- `Phosphorylation_site_dataset.gz`

The file version marker is `042026`. It contains `ACC_ID`, `ORGANISM`, and
`MOD_RSD`, which are sufficient to compute human phosphorylation-site counts.
Preserve PhosphoSitePlus attribution and license requirements in any public
release or derived resource.
