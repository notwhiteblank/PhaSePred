# IUPred3

Local archive currently present:

- `iupred3.tar.gz`

This is a local IUPred3 source release for protein disorder prediction. The
archive and extracted tool tree are ignored by git; this README records the
local setup.

## Local Install

Extract in place:

```bash
tar -xzf tools/per-tool/IUPred3/iupred3.tar.gz -C tools/per-tool/IUPred3
```

Expected extracted tree:

- `tools/per-tool/IUPred3/iupred3/iupred3.py`
- `tools/per-tool/IUPred3/iupred3/iupred3_lib.py`
- `tools/per-tool/IUPred3/iupred3/data/`
- `tools/per-tool/IUPred3/iupred3/P53_HUMAN.seq`

## Runtime

IUPred3 is a Python 3 command-line tool and requires SciPy. In this workspace it
runs with the project `uv` environment.

Use the project wrapper from the repository root:

```bash
tools/wrappers/run_iupred3.sh input.fasta long
```

Prediction types:

- `long`: long disorder
- `short`: short disorder
- `glob`: globular domain prediction

Optional flags:

- `-a` or `--anchor`: include ANCHOR2 prediction
- `-s no|medium|strong`: smoothing mode; upstream default is `medium`

## Smoke Test

```bash
tools/wrappers/run_iupred3.sh tools/per-tool/IUPred3/iupred3/P53_HUMAN.seq long \
  | sed -n '1,20p'
```

Expected output includes the IUPred3 citation header, the selected prediction
type, and a per-residue table headed `# POS RES IUPRED2`.
