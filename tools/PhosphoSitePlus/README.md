# PhosphoSitePlus data contract

PhaSePred's human modes (hSaPS / hPdPS) use the PhosphoSitePlus
phosphorylation-site dataset for the **Phos freq** feature column. This contract
locates, installs, and verifies that dataset; it is a data contract
(`kind = "data"`), not a runnable tool.

## Install

```bash
bash tools/PhosphoSitePlus/install.sh            # download + verify
bash tools/PhosphoSitePlus/install.sh --check    # verify only (exit 0/1)
bash tools/PhosphoSitePlus/install.sh --offline  # never use the network
```

The installer writes only into the user data directory:

| Tier | Location |
|---|---|
| `$PHASEPRED_DATA_ROOT` set | `$PHASEPRED_DATA_ROOT/phosphositeplus/` |
| else `$XDG_DATA_HOME` set | `$XDG_DATA_HOME/phasepred/phosphositeplus/` |
| else | `~/.local/share/phasepred/phosphositeplus/` |

A legacy checkout copy at `data/raw/external/phosphositeplus/` is still found by
the resolver (tier 3) and reported as `source = "repo"`.

To install from a pre-fetched file instead of the network, place
`Phosphorylation_site_dataset.gz` in `$PHASEPRED_VENDOR_ARCHIVE_DIR`; the
installer copies it and still verifies the digest.

## Integrity policy — known-version, not hard

PhosphoSitePlus is a **rolling database**: upstream revises the file without a
version bump, so a pinned sha256 cannot be a hard gate or the installer would
break within months. The manifest therefore declares:

```toml
sha256 = "5e339a7f…b095"          # measured 2026-09-17, version marker 042026
sha256_policy = "known-version"
```

- **match** → silent success.
- **mismatch** → a warning naming both digests, the actual digest is written to
  `phosphositeplus/install-record.txt`, and installation still exits 0.

This is the deliberate opposite of **PScore**, whose 2019 static archive uses
`sha256_policy = "hard"` and aborts on mismatch. See `docs/TOOL_LICENSES.md` and
the E4 install documentation.

## Licence

**CC BY-NC-SA 3.0 Unported** (Creative Commons
Attribution-NonCommercial-ShareAlike 3.0). The dataset is
**not redistributable** — it is fetched at install time and never committed. If
you use it, the dataset's own notice requires:

- the words "PhosphoSitePlus(R), www.phosphosite.org" in the text or webpage, and
- the citation: Hornbeck PV, Zhang B, Murray B, Kornhauser JM, Latham V,
  Skrzypek E. *PhosphoSitePlus, 2014: mutations, PTMs and recalibrations.*
  Nucleic Acids Res. 2015;43:D512-20. PMID: 25514926.

## Usage from phasepred

`phasepred check-tools` reports the row as `PhosphoSitePlus` (`source = "user"`
for a user-root copy, `"repo"` for a checkout copy). `phasepred predict` uses it
through `phasepred.tools.phosphosite_path()`; when it is absent, human-mode
prediction proceeds with median imputation and emits
`PhaSePredMissingDataWarning` rather than aborting.
