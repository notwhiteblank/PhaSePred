# ESpritz — tool package

Self-contained ESpritz contract (PhaSePred S3 tool standardization).

## What is ESpritz

ESpritz (Tosatto lab) predicts disordered residues. PhaSePred uses the
**DisProt** model `D` at the **5% false-positive-rate threshold** (`sw 0`)
for the `IDR` feature (fraction of `D`-state residues; banner/license lines
are ignored).

- Upstream: <https://protein.bio.unipd.it/espritz/>
- License: **Tosatto lab academic license v1.1 — non-commercial,
  no redistribution** (`tools/ESpritz/espritz/LICENSE`). See
  `docs/TOOL_LICENSES.md` §5.
- Reference: Walsh I., Martin A.J.M., Di Domenico T., Tosatto S.C. (2012)
  ESpritz: accurate and fast prediction of protein disorder.
  *Bioinformatics* 28:503-509.

## Package layout

```
tools/ESpritz/
├── run            entrypoint: ./run <fasta_dir> <MODEL> <SW>
├── manifest.toml  contract (kind=perl, license=Tosatto-academic, NOT redistributable)
├── install.sh     downloads-guide + extract (interactive; --check for verify-only)
├── README.md      this file
├── espritz.zip    gitignored user-supplied archive
└── espritz/       gitignored extracted tree
```

ESpritz cannot be redistributed: only the contract three files are tracked;
the archive and extracted tree stay out of git.

## Install

```bash
bash tools/ESpritz/install.sh          # interactive: waits for espritz.zip
bash tools/ESpritz/install.sh --check  # verify already-extracted tree only
```

## Usage

```bash
tmpdir="$(mktemp -d)"; cp in.fasta "$tmpdir/"
tools/ESpritz/run "$tmpdir" D 0     # paper protocol (DisProt, 5% FPR)
tools/ESpritz/run --check
```

`MODEL` ∈ {X, D, N} (sequence-only models; pX/pD/pN need PSI-BLAST and are
blocked). `SW` ∈ {0 = 5% FPR, 1 = best-Sw}. Writes `<acc>.espritz` sibling
files with two-column residue lines (declared output format in
`manifest.toml`: `espritz`).

## Attribution / license notice

The run banner shows a licensed-to notice; do not remove it:

> Licensed to: PhD Kaiqiang You (Academic, Biomedical infomatics, Peking
> University) — this license is for non-commercial use only.

Keep `analysis/tools/ESpritz/espritz/LICENSE` intact whenever the extracted
tree is present.