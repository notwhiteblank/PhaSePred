# PhaSePred training tables (`data/processed/`)

This directory holds the four training tables for the PhaSePred 2022-protocol
baseline (E1). The format is owned exclusively by
`src/phasepred/training_data.py`; the build driver is
`scripts/build_training_tsv.py`. Nothing here is hand-edited.

## Provenance

The tables are a faithful serialization of Chen et al. 2022 *PNAS*
`pnas.2115369119.sd02.xlsx` (training) and `pnas.2115369119.sd03.xlsx`
(independent test), read through `phasepred.paper.load_paper_features()`
(header row 2, `header=1` fallback). The two workbooks are **not
redistributed**. Their measured sha256 digests are
recorded in `MANIFEST.tsv`:

- `pnas.2115369119.sd02.xlsx` — `47675052810146ba62a60910ed37db52475dd6d0bae64b3ee3066fa417631cad`
- `pnas.2115369119.sd03.xlsx` — `298f5293397298f3647c5ddeeb104edee243ed23752883e2bf653a0b40e69ddf`

## Files (scope × split)

Measured values transcribed from `MANIFEST.tsv`:

| file | scope | split | rows | columns | bytes | sha256 |
|---|---|---|---|---|---|---|
| `chen2022_s2s3_base-features_train.tsv` | base | train | 96658 | 17 | 13399739 | `df9081b001936adb678de36ac186811df7ed0f5db1e09079e6e788049d385bb7` |
| `chen2022_s2s3_base-features_test.tsv` | base | test | 24416 | 18 | 3508476 | `b26253f921217375f4873d5714058fcf20925729b21b8381a4eba27f3940ff9a` |
| `chen2022_s2s3_human-features_train.tsv` | human | train | 17757 | 19 | 2800875 | `fedecb3850263de9b561a6a5685f91b3f1235135b4baf6aca4504031969cad59` |
| `chen2022_s2s3_human-features_test.tsv` | human | test | 4540 | 20 | 736518 | `cb36429d16f575785213dc625708c302e093129e4dc29ddad60facad36928134` |

`base` covers the non-human models (`SaPS`, `PdPS`); `human` covers the human
models (`hSaPS`, `hPdPS`). Test tables carry the leading `Source` column;
train tables do not.

## Column semantics

Shared columns:

| column | meaning |
|---|---|
| `UniprotEntry` | UniProt accession (string). |
| `task` | `SaPS` / `PdPS` (base) or `hSaPS` / `hPdPS` (human). |
| `split` | `train` or `test`. |
| `label` | `1` = positive, `0` = negative. |
| `source_sheet` | workbook sheet the row came from (`SaPS`, `NoPS`, `SaPS-test`, …). |
| `Source` | **test tables only**; provenance database (`PhaSepDB`, `LLPSDB`, `PhaSePro`); empty for negatives, whose source sheets lack the column. |
| `Gene name` | gene symbol from the paper sheet. |
| `Organism` | source organism. |
| `Organism ID` | NCBI taxonomy id. |
| `length` | peptide length in residues. |

Feature columns, in schema order:

| column | scope | definition |
|---|---|---|
| `Hydropathy` | both | Kyte–Doolittle hydropathy. |
| `FCR` | both | fraction of charged residues (D/E/K/R); §2. |
| `IDR` | both | ESpritz-DisProt disordered fraction; §3. |
| `LCR` | both | SEG low-complexity fraction; §4. |
| `PScore` | both | PScore phase-separation score; §5. |
| `PLAAC` | both | normalized PLAAC LLR; §6. |
| `catGRANULE` | both | catGRANULE score; §7. |
| `DeepCoil` | both | binarized coiled-coil indicator (0/1, threshold 0.82); §8. |
| `Phos freq` | human | PhosphoSitePlus sites / length; §9. |
| `DeepPhase` | human | DeepPhase droplet-forming probability; §10. |

The authoritative feature definitions are not repeated here.

## Missing values

A missing value is written as an **empty field**; there is no sentinel
(`NaN`, `<NA>`, `NULL`, `None`, `nan` never appear), and **no imputation is
applied**. This is deliberate: the paper protocol passes `NaN` straight to
XGBoost's built-in missing-value handling. To read back, use
`phasepred.training_data.read_table()`, whose `na_values=[""]` +
`float_precision="round_trip"` preserve both missingness and 17-digit float
precision exactly.

Measured per-column missing counts (`DataFrame.isna().sum()`):

### base / train (96658 rows)

```
UniprotEntry     0    Gene name    5    Hydropathy  160
task             0    Organism     0    FCR           0
split            0    Organism ID  0    IDR           0
label            0    length       0    LCR           0
source_sheet     0                       PScore     9302
                                        PLAAC       412
                                        catGRANULE    7
                                        DeepCoil      0
```

### base / test (24416 rows)

```
Source       24124    Gene name    1    Hydropathy   45
UniprotEntry     0    Organism     0    FCR           0
task             0    Organism ID  0    IDR           0
split            0    length       0    LCR           0
label            0                       PScore     2260
source_sheet     0                       PLAAC        96
                                        catGRANULE    23
                                        DeepCoil       2
```

### human / train (17757 rows)

```
UniprotEntry     0    Gene name    0    Hydropathy   30
task             0    Organism     0    FCR           0
split            0    Organism ID  0    IDR           0
label            0    length       0    LCR           0
source_sheet     0                       PScore     1662
                                        PLAAC        60
                                        catGRANULE    4
                                        DeepCoil      0
                                        Phos freq     0
                                        DeepPhase  7222
```

### human / test (4540 rows)

```
Source        4400    Gene name    0    Hydropathy    4
UniprotEntry     0    Organism     0    FCR           0
task             0    Organism ID  0    IDR           0
split            0    length       0    LCR           0
label            0                       PScore      392
source_sheet     0                       PLAAC        12
                                        catGRANULE    4
                                        DeepCoil      0
                                        Phos freq     0
                                        DeepPhase  1874
```

## Duplicate rows

The dedicated negative sheets (`NoPS`, `hNoPS`) are shared by two tasks each,
so the long-form tables contain cross-task duplicates of the same accession
with different `task` values. These are part of the protocol and are preserved.
There are **no duplicates within a single `(task, split)`**, which
`validate_schema` enforces.

Measured unique-accession counts and duplication (`rows − nunique(UniprotEntry)`):

| table | rows | unique accessions | duplicate rows |
|---|---|---|---|
| base / train | 96658 | 48500 | 48158 |
| base / test | 24416 | 12301 | 12115 |
| human / train | 17757 | 8956 | 8801 |
| human / test | 4540 | 2317 | 2223 |

## Row order

Row order and row multiplicity are copied exactly from
`load_paper_features()`: no sorting and no de-duplication. Training samples by
positional index (`numpy.random.default_rng(42 + i)`), so row order is part of
the protocol and must not be rearranged or regenerated with a different
`load_paper_features()`.

## Reproduction

```
python scripts/build_training_tsv.py            # build (needs the private workbooks; takes minutes)
python scripts/build_training_tsv.py --check    # verify existing tables (no workbooks needed)
```

The build is deterministic: no wall-clock value enters any output, and two
consecutive builds are byte-identical. `SCRIPT_VERSION` is bumped whenever a
serialization change alters the output bytes. The audited sizes live in
`EXPECTED_ROWS` / `EXPECTED_BYTES`, and the frozen digests in `EXPECTED_SHA256`;
`validate_expected_size()` fails on any byte, digest or row-count drift, so the
constants act as a drift detector. Changing the tables deliberately requires:
bump `SCRIPT_VERSION`, regenerate, re-freeze `EXPECTED_SHA256`, update
`EXPECTED_BYTES`.

## Schema enforcement

`phasepred.training_data` is the single authority for the on-disk format. A
`TrainingDataSchemaError` is raised when:

1. an unknown `scope` is requested (`_key`);
2. an unknown `split` is requested (`_key`);
3. the column list or column order differs from `COLUMNS[(scope, split)]`
   (`validate_schema`);
4. the file contains a task outside `TASKS_BY_SCOPE[scope]` (`validate_schema`);
5. the `split` column is not the single expected split (`validate_schema`);
6. a `label` is not `0` or `1` (`validate_schema`);
7. a `(task, split, UniprotEntry)` triple is duplicated (`validate_schema`);
8. a source frame lacks a required schema column (`select_and_cast`);
9. a table file is missing (`read_table`, `validate_expected_size`);
10. the byte count differs from `EXPECTED_BYTES` (`validate_expected_size`);
11. the digest differs from `EXPECTED_SHA256` once frozen
    (`validate_expected_size`);
12. the row count differs from `EXPECTED_ROWS` (`validate_expected_size`);
13. `MANIFEST.tsv` is missing (`read_manifest`);
14. a task is requested that does not belong to the scope (`load_task_frame`).
