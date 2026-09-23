# External Tools

PhaSePred's prediction pipeline calls into six external feature tools, plus
one data-only contract (`tools/PhosphoSitePlus/`) that pins a download rather
than a runnable tool. Since S3 every tool is a self-contained **tool package**
under `tools/<TOOL>/` with a machine-readable contract:

- `run` — unified entrypoint (`./run <in.fasta> [native args]` → stdout in
  the `manifest.toml`-declared `output` format; `./run --check` verifies
  the runtime entity in place).
- `manifest.toml` — `name/kind/entry/license/redistributable/runtime/
  path_executables/data/output/install` (schema in `docs/FEATURES.md`
  附注 A and `docs/TOOL_LICENSES.md` §12).
- `install.sh` — idempotent, non-interactive installer; `--check` verifies
  only and `--offline` forbids the network. Downloads land in the user data
  root, never in the checkout.
- `README.md` — usage / provenance / license / attribution.

phasepred locates tools through the manifest registry in
`src/phasepred/tool_paths.py` (`Runner` reads the contracts; it no longer
hardwires any per-tool path).

## Tool packages

| Package | kind | License (SPDX) | Redistributable | Runtime |
|---|---|---|---|---|
| `tools/SEG/` | binary | NCBI public-domain convention | yes | local build (`make`) |
| `tools/PLAAC/` | java | MIT | yes | `java>=11` |
| `tools/PScore/` | python | CC-BY-4.0 | yes (code); DBS data gitignored | `python>=3.12` + numpy |
| `tools/ESpritz/` | perl | Tosatto academic v1.1 | **no** | `perl` |
| `tools/DeepCoil/` | env-prefix | none declared | **no** | isolated conda env |
| `tools/LocalCIDER/` | python-library | GPL-2.0 | yes (via pip, not vendored) | project Python |
| `tools/PhosphoSitePlus/` | data | CC-BY-NC-SA-3.0 | **no** | none (registration-gated dataset) |

`tools/_archive/` holds historical tools: `catGRANULE_v2/` (53/82 scales
defect — not in the predict path) and `IUPred3/` (paper does not use it,
D17). Neither is wired into phasepred.

## Resolution cascade

For every tool the manifest resolver in `src/phasepred/tool_paths.py`
tries, in order:

1. **Environment variable** `PHASEPRED_<TOOL>_DIR` (e.g.
   `PHASEPRED_SEG_DIR`, `PHASEPRED_ESPRITZ_DIR`).
2. **User data root / repo checkout** — for contracts registered in
   `data._COMPONENT_LEGACY_NAMES` (PScore, ESpritz, PhosphoSitePlus), a
   complete package under `$PHASEPRED_DATA_ROOT` / `$XDG_DATA_HOME/phasepred`
   / `~/.local/share/phasepred`, or the checkout's own
   `data/raw/external/<name>/`.
3. **Vendored tool package** `tools/<TOOL>/` in this repository.
4. **PATH** (`shutil.which`) — only where the manifest `path_executables`
   lists a sensible name (SEG → `seg`, DeepCoil → `deepcoil`).
5. Error with the `install.sh` hint from the manifest.

`check-tools` reports which tier supplied each hit through `source`, using the
same vocabulary as `docs/INSTALL.md`: `env`, `user` (a writable user root),
`repo` (the checkout fallback), `vendored`, `path`, `package`, `missing`.
DeepCoil (`kind = "env-prefix"`) is the exception: its isolated env resolves
through `$DEEPCOIL_ENV_PREFIX` → user root `envs/deepcoil` → legacy
`.external_envs/deepcoil`. See `docs/INSTALL.md` §1 and §4 for the full cascade
and its tier definitions.

`check-tools` additionally runs each package's `./run --check` so the
report reflects *runtime usability* (binary built / data extracted / env
installed), not just the presence of the contract files.

Run `uv run phasepred check-tools` after installing to verify.

## Install overview

```bash
bash tools/SEG/install.sh                  # builds seg (or make in tools/SEG)
bash tools/PScore/install.sh               # downloads + hard-verifies + extracts
bash tools/PLAAC/install.sh                # no-op when jar present
bash tools/ESpritz/install.sh              # needs a vendor archive to be non-interactive
bash tools/DeepCoil/install.sh             # isolated conda env + pip
bash tools/LocalCIDER/install.sh           # pip install localcider
bash tools/PhosphoSitePlus/install.sh      # hSaPS/hPdPS data lookup
```

Each supports `--check` (verify only, exit 0/1) and `--offline` (skip
downloads; exit 0 when already satisfied, non-zero otherwise). `phasepred
predict` aborts up-front if anything required is MISSING; the full install
matrix is in `docs/INSTALL.md`.