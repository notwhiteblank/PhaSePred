# PScore

Local archive currently present:

- `elife-31486-code2-v2.tgz`

This is the eLife supporting archive `Source code 2` for the PScore feature
tool. The archive and extracted source tree are ignored by git.

## Local Install

Extract in place:

```bash
tar -xzf tools/per-tool/PScore/elife-31486-code2-v2.tgz -C tools/per-tool/PScore
chmod +x tools/per-tool/PScore/SourceCodeS2/elife_phase_separation_predictor.py
```

Expected extracted tree:

- `tools/per-tool/PScore/SourceCodeS2/elife_phase_separation_predictor.py`
- `tools/per-tool/PScore/SourceCodeS2/DBS/`
- `tools/per-tool/PScore/SourceCodeS2/PDB.TESTSET.fasta`

## Runtime

The script requires Python plus `numpy`. It was smoke-tested in this repository
with:

- Python 3.12.12 through `uv run python`
- numpy 2.4.4

The original script expects to be run from the directory containing `DBS/`,
unless its internal `DBPATH` variable is changed.

## Usage

Run from the extracted source directory:

```bash
cd tools/per-tool/PScore/SourceCodeS2
uv run python elife_phase_separation_predictor.py input.fasta \
  -output pscore.tsv \
  -overwrite
```

Useful options from the bundled README:

- `-residue_scores`: write weighted per-residue scores.
- `-score_components`: write the eight unweighted per-residue score components.
- `-output FILE`: write to a file instead of stdout.
- `-overwrite`: replace an existing output file.
- `-mute`: suppress progress messages in output-file mode.

## Input And Output

Input is a FASTA file. The tool skips sequences shorter than the supported
training/test-set range and sequences containing ambiguous residues.

Default output is one line per scored sequence:

```text
PScore:                2.03                  >|1a12A
```

## Smoke Test

The bundled `PDB.TESTSET.fasta` contains thousands of sequences and is slow as a
full test. For a quick end-to-end check, run PScore against the first bundled
FASTA record and write output outside the repo:

```bash
cd tools/per-tool/PScore/SourceCodeS2
uv run python elife_phase_separation_predictor.py \
  <(awk 'BEGIN{n=0} /^>/{n++; if(n>1) exit} {print}' PDB.TESTSET.fasta) \
  -output /tmp/pscore_quick.tsv \
  -overwrite \
  -mute
sed -n '1,8p' /tmp/pscore_quick.tsv
```

Expected smoke-test output begins:

```text
PScore:                2.03                  >|1a12A
```
