# PLAAC

Local archive currently present:

- `plaac-master.zip`

This is a local GitHub source archive for PLAAC. The archive and extracted
source tree are ignored by git.

## Local Install

Extract in place:

```bash
unzip tools/per-tool/PLAAC/plaac-master.zip -d tools/per-tool/PLAAC
```

Expected extracted tree:

- `tools/per-tool/PLAAC/plaac-master/cli/`
- `tools/per-tool/PLAAC/plaac-master/web/bin/plaac.jar`
- `tools/per-tool/PLAAC/plaac-master/LICENSE.TXT`

## Runtime

The command-line implementation is Java. This machine currently has a Java
runtime available, but not the JDK build tools `javac` and `jar`, so the
precompiled upstream jar at `web/bin/plaac.jar` is the usable local CLI artifact.

`Rscript` is available locally for the optional plotting scripts.

## Usage

Use the project wrapper from the repository root:

```bash
tools/wrappers/run_plaac.sh -i input.fasta > plaac.tsv
```

For per-residue plot data:

```bash
tools/wrappers/run_plaac.sh -i input.fasta -p all > plaac_plotdata.tsv
```

PLAAC defaults to Saccharomyces cerevisiae background frequencies
(`-a 1.0`). For non-yeast or small candidate FASTA files, provide a proteome
background with `-b background.fasta` or a frequency table with `-B bg_freqs.txt`
before using interpolated backgrounds such as `-a 0.5`.

## Optional Build

If a JDK is installed, rebuild the CLI jar from source with:

```bash
cd tools/per-tool/PLAAC/plaac-master/cli
./build_plaac.sh
```

The expected built jar is `tools/per-tool/PLAAC/plaac-master/cli/target/plaac.jar`.

## Smoke Test

```bash
tools/wrappers/run_plaac.sh \
  -i tools/per-tool/PLAAC/plaac-master/cli/example/four_classic_prions.fasta \
  | sed -n '1,20p'
```

Expected output includes the PLAAC runtime parameter block followed by a tabular
header beginning with `SEQid`.
