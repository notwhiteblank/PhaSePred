# Changelog

## 1.0.1 (2026-09-16)

- **Version declarations converge.** `pyproject.toml` and
  `catgranule.__version__` now both say 1.0.1 (E2, D38).
- **PyPI correction.** The 1.0.0 entry's claim that the package was "not
  uploaded to PyPI" was wrong: catgranule 1.0.0 was uploaded on
  2026-09-07 (see the internal migration plan §9, S10 errata).
- **xgboost constraint unchanged** (`>=3.2,<4`). The distilled `.ubj`
  Booster yields bit-identical `single` and `granule_strength` values under
  xgboost 3.2.0 and 3.4.1 (first measured 2026-09-16; re-measured on
  Q08211/P35637/O75146/O95153 on 2026-09-23, E7), so no upper-bound
  tightening is warranted. The `.ubj` Booster is shipped and never retrained
  by this project, so only inference runs, and inference is insensitive to
  the xgboost minor version.
- **Upload status.** As of 2026-09-23, catgranule 1.0.1 has not been uploaded
  to PyPI; 1.0.0 (2026-09-07) remains the latest published version. The
  intended upload target is **1.0.1**, before `phasepred` 1.0.0 (E9).

## 1.0.0 (2026-09-07)

- **Default weights: distilled surrogate (route=distilled).** XGBoost model
  reproduced from the locked S4 configuration (600 trees / depth 7 /
  lr 0.04, random_state 0, 181-dim feature map, trained on the web-archive
  human eval-pool only), accepted by user decision D21 (gap-acceptance;
  holdout `rho`=0.998961, median `|delta|`=0.0177). The 1.0.0 reproduction
  is gate-exact (see `docs/CATGRANULE_VALIDATION.md` S4b section and
  `scripts/reproduce_distilled.py`).
- **Dual-provenance weights.** Shipped artifacts now include:
  - `catgranule_distilled_weights.json` + `distilled_xgb_model.ubj`
    (default; `mode="distilled"`, self-contained feature basis).
  - `catgranule_paper_weights.json` (alternate; Bolognesi 2016 closed-form
    recipe with S4-recomputed yeast Z constants; rho 0.9747 / delta 0.186).
  - `catgranule_v1_weights.json` (legacy-unaudited; retained for
    continuity).
- **Weight switching.** `CATGRANULE_WEIGHTS` env var, `route=` API
  parameter, and the CLI `--route` / `--weights` options select the
  artifact. Distilled is the default; both paper and legacy remain
  switchable and fully tested.
- **API compatibility.** `score_sequence` / `score_batch` signatures
  unchanged (added optional `route`). CLI `predict` subcommand outputs
  `single` + `residue` + `granule_strength`; new `check` output includes the
  active artifact mode.
- **Residue profile.** The 51-window / offset-24 / L-50 geometry is
  preserved on all routes; the distilled route aligns the profile to the
  surrogate `single` (mean-equal per protein).
- **Dependencies.** `xgboost>=3.2,<4`, `numpy>=1.26,<3` and `typer>=0.15`
  declared (runtime scoring uses the pure `xgboost.Booster` API — no
  scikit-learn dependency; no external tools).
- Package uploaded to PyPI as `catgranule` 1.0.0 on 2026-09-07 (D18
  settled; internal migration plan §9, S10 errata).
