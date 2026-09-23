# PLAAC — tool package

Self-contained PLAAC contract (PhaSePred S3 tool standardization).

## What is PLAAC

PLAAC (Prion-Like Amino Acid Composition, Whitehead Institute) predicts
prion-like amyloid-forming composition. PhaSePred uses the **NLLR**
(normalized log-likelihood ratio) column for the `PLAAC` 8-feature model.

- Upstream: <http://plaac.wi.mit.edu/> (GitHub `whitehead/plaac`)
- License: **MIT** (see `plaac-master/LICENSE.TXT`) — redistributable.
- Reference: Lancaster, A.K., Nutter-Upham, A., Lindquist, S., King, O.D.
  (2014) PLAAC: a web and command-line application to identify proteins
  with prion-like amino acid composition. *Bioinformatics* 30:2501-2502.

## Package layout

```
tools/PLAAC/
├── run            entrypoint: ./run -i in.fasta [other plaac args]
├── manifest.toml  contract (kind=java, license=MIT, redistributable)
├── install.sh     verifies/extracts the vendored tree (--check supported)
├── README.md      this file
└── plaac-master/  vendored upstream tree (LICENSE.TXT + cli/src/ + web/bin/plaac.jar)
```

`plaac-master.zip` (large upstream archive) is gitignored; the extracted
source and the prebuilt `web/bin/plaac.jar` are tracked.

## Install

```bash
bash tools/PLAAC/install.sh          # no-op when jar present
bash tools/PLAAC/install.sh --check  # verify only
```

Requires a Java runtime (`openjdk >= 11`; 11.0.32 verified on the dev host,
see `docs/TOOL_LICENSES.md` §11).

## Usage

```bash
tools/PLAAC/run -i input.fasta > plaac.tsv
tools/PLAAC/run --check
```

Standard output: PLAAC runtime parameter block + tab-separated table
(header `SEQid` → `.tsv`). PhaSePred reads the `NLLR` column
(declared output format in `manifest.toml`: `tsv`).

PLAAC defaults to *Saccharomyces cerevisiae* background frequencies
(`-a 1.0`). For non-yeast or small candidate FASTA files, provide a
proteome background with `-b background.fasta` or a frequency table with
`-B bg_freqs.txt` (see the vendored `cli/README.md`).

## Rebuilding the jar (optional, needs JDK)

```bash
cd tools/PLAAC/plaac-master/cli && ./build_plaac.sh
```

## Attribution

MIT license: copyright (c) 2009-2014 Whitehead Institute for Biomedical
Research, 2011 Boston Biomedical Research Institute, 2013 University of
Massachusetts Medical School. Keep `plaac-master/LICENSE.TXT` intact.