# SEG — tool package

Self-contained SEG contract (PhaSePred S3 tool standardization).

## What is SEG

SEG is the NCBI low-complexity masking program (Wootton & Federhen 1993).
The PhaSePred paper uses it for the **LCR** feature: the fraction of
residues masked as low-complexity under the default parameters
(`12 2.2 2.5 -x`).

- Source: <https://ftp.ncbi.nih.gov/pub/seg/seg/> (downloaded 2026-05-07)
- License: no LICENSE file — US-government / public-domain convention
  (see `docs/TOOL_LICENSES.md` §2). Redistributable.
- Original `README` + `seg.doc`, the source (`*.c`, `*.h`, `*fac.h`) and
  the `makefile` are vendored in this directory and tracked by git.

## Package layout

```
tools/SEG/
├── run            entrypoint: ./run <in.fasta> [12 2.2 2.5 -x]
├── manifest.toml  machine-readable contract (name/kind/license/install/...)
├── install.sh     idempotent build (--check verifies the binary)
└── README.md      this file
```

Entities for the binary built from source: `seg.c`, `genwin.c`, `genwin.h`,
`hiseg.c`, `lnfac.h`, `makefile`, `README`, `seg.doc`, plus the prebuilt
Linux `seg` binary. `seg.sgi.iris5`/`seg.sun.sunos4` and `*.o` are
gitignored local build artifacts.

## Install

```bash
bash tools/SEG/install.sh          # builds ./seg (idempotent)
bash tools/SEG/install.sh --check  # verify only
```

## Usage

```bash
tools/SEG/run tests/fixtures/sequences/Q08211.fasta            # paper args applied
tools/SEG/run in.fasta 12 2.2 2.5 -x                          # explicit args
tools/SEG/run --check                                         # status probe
```

Standard output = masked FASTA (headers preserved, masked residues as
lowercase / `x` under `-x`). Declared output format in `manifest.toml`:
`text`.

## Attribution

Keep the NCBI `README` (with the Wootton & Federhen 1993 Computers and
Chemistry citation) next to the source. Reference:

> Wootton, J. C. and S. Federhen (1993). Statistics of local complexity in
> amino acid sequences and sequence databases. *Computers and Chemistry*
> 17:149-163.