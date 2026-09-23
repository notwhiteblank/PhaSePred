# PScore — tool package

Self-contained PScore contract (PhaSePred S3 tool standardization).

## What is PScore

PScore is the pi-pi-contact residue-interaction predictor from the Forman-Kay
lab. PhaSePred uses its **single PScore value** for the `PScore` 8-feature
model.

- Upstream: eLife article 31486, "Source code 2"
  (`elife-31486-code2-v2.tgz`)
- License: **CC-BY 4.0** (eLife article and supplementary materials;
  program files carry no separate license header) — redistributable with
  attribution. See `docs/TOOL_LICENSES.md` §4.
- Reference: Vernon, R.M., Chong, P.A., Tsang, B., Kim, T.H., Bah, A.,
  Farber, P., Lin, H., Forman-Kay, J.D. (2018) Pi-Pi contacts are an
  overlooked protein feature relevant to phase separation.
  *eLife* 7:e31486. DOI 10.7554/eLife.31486

## Package layout

```
tools/PScore/
├── run            entrypoint: ./run <in.fasta> [-output ... -overwrite -mute]
├── manifest.toml  contract (kind=python, license=CC-BY-4.0)
├── install.sh     extracts SourceCodeS2/ from the eLife archive (--check)
├── README.md      this file
├── SourceCodeS2/      vendored code (tracked): predictor + README.txt + PDB.TESTSET.fasta
│   └── DBS/           gitignored data (135 MB): grids, freq tables, WeightSettings.txt
└── elife-31486-code2-v2.tgz   gitignored upstream archive
```

The **code** is vendored and tracked; the **DBS/ data** (135 MB) is kept
gitignored (deferred upgrade; can be re-elevated when desired, D15).

## Install

```bash
bash tools/PScore/install.sh          # extract SourceCodeS2/ + DBS/ (idempotent)
bash tools/PScore/install.sh --check  # verify predictor + DBS present
```

## Usage

Run from anywhere; the wrapper `cd`s into `SourceCodeS2/` because the
predictor reads `DBS/` relative to its own location:

```bash
tools/PScore/run input.fasta -output pscore.tsv -overwrite -mute
tools/PScore/run --check
```

Useful options (from the vendored README.txt): `-residue_scores`,
`-score_components`, `-output FILE`, `-overwrite`, `-mute`. Sequences
shorter than the supported range or with ambiguous residues are skipped.

Standard output: one `PScore:  <value>  >ACC` line per scored accession
(declared output format in `manifest.toml`: `text`).

## Attribution (CC-BY 4.0)

Carry the reference above and the CC-BY 4.0 notice. The archived
`SourceCodeS2/README.txt` documents the original authors. Data derived from
the upstream `DBS/` must be attributed to Vernon et al. (2018).