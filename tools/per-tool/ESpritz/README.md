# ESpritz

Local archive currently present:

- `espritz.zip`

This is a local ESpritz 1.1 package for protein disorder prediction. The archive
and extracted tool tree are ignored by git; this README records the local setup.

## Local Install

Extract in place:

```bash
unzip tools/per-tool/ESpritz/espritz.zip -d tools/per-tool/ESpritz
chmod +x tools/per-tool/ESpritz/espritz/espritz.pl tools/per-tool/ESpritz/espritz/bin/disbin*
chmod +x tools/per-tool/ESpritz/espritz/bin/blastredo
```

Expected extracted tree:

- `tools/per-tool/ESpritz/espritz/espritz.pl`
- `tools/per-tool/ESpritz/espritz/bin/disbinD`
- `tools/per-tool/ESpritz/espritz/bin/disbinN`
- `tools/per-tool/ESpritz/espritz/bin/disbinX`
- `tools/per-tool/ESpritz/espritz/models/`
- `tools/per-tool/ESpritz/espritz/example_fastas/`

## Runtime

The sequence-only modes use Perl plus bundled Linux x86-64 binaries. ESpritz
expects to run from its own installation directory, so use the project wrapper
from the repository root:

```bash
tools/wrappers/run_espritz.sh input_fasta_dir D 0
```

Arguments:

- `MODEL`: `X` for X-ray disorder, `D` for DisProt disorder, or `N` for NMR
  disorder.
- `SW`: `1` for the best-Sw threshold, or `0` for the 5% false-positive-rate
  threshold.

Input files must have the `.fasta` extension. ESpritz writes sibling output
files ending in `.espritz` and `.espritz.fasta`, and it rewrites temporary
normalized FASTA files while running.

For PhaSePred paper-style IDR feature reconstruction, use ESpritz DisProt at the
5% false-positive-rate threshold:

```bash
tools/wrappers/run_espritz.sh input_fasta_dir D 0
```

The PhaSePred article describes "the ESpritz DisProt program" with the decision
threshold set at 5% false-positive rate, and does not mention PSI-BLAST/PSSM
input for this feature.

## PSI-BLAST Modes

The upstream `pX`, `pD`, and `pN` modes require additional BLAST binaries,
sequence databases, and edits to `espritz/align/getAlignments.pl`. The project
wrapper intentionally blocks those modes until the external database/toolchain
layer is documented.

## Smoke Test

Use a copy of the bundled examples so the ignored upstream example files are not
modified:

```bash
tmpdir="$(mktemp -d)"
cp tools/per-tool/ESpritz/espritz/example_fastas/*.fasta "$tmpdir/"
tools/wrappers/run_espritz.sh "$tmpdir" D 0
ls "$tmpdir"/*.espritz "$tmpdir"/*.espritz.fasta
```

Expected output includes ESpritz's license banner and one `.espritz` plus one
`.espritz.fasta` file per input `.fasta`.
