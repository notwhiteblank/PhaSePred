# SEG

Local mirror of the NCBI SEG low-complexity masking tool.

Source URL: <https://ftp.ncbi.nih.gov/pub/seg/seg/>

Downloaded on 2026-05-07T08:43:50Z.

## Local Files

- `README`
- `genwin.c`
- `genwin.h`
- `hiseg.c`
- `lnfac.h`
- `makefile`
- `seg.c`
- `seg.doc`
- `seg.sgi.iris5`
- `seg.sun.sunos4`

The source files, downloaded binaries, local build outputs, and this directory
tree are ignored by git except for this README.

## Build

```bash
cd tools/per-tool/SEG
make
```

This builds `./seg` from the NCBI source. The code is old K&R-style C and
modern compilers emit many implicit declaration warnings, but the local build
completed with the bundled makefile on 2026-05-07.

## Usage

```bash
tools/per-tool/SEG/seg <input.fasta> 12 2.2 2.5 -x
```

Run `tools/per-tool/SEG/seg` with no arguments to show the bundled usage text. See
`seg.doc` for the original manual page.
