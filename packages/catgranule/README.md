# catgranule

catGRANULE phase-separation propensity scoring, following Bolognesi et al.
2016, *Cell Reports* 16:222-231 ("A concentration-dependent liquid phase
separation can cause toxicity upon increased protein expression").

This package re-implements the paper's scoring equations (Equation 1-2:
six normalized physico-chemical propensities averaged over a centered
heptapeptide window, plus a log-length term; Equation 3: granule
strength). It is the scoring engine used by the
[PhaSePred](https://www.pnas.org/doi/10.1073/pnas.2115369119) re-build
(phasepred package).

## Scoring

```python
from catgranule import score_sequence, score_batch

result = score_sequence("MGGYNNNNSS...")          # default: distilled surrogate
# {"single": ..., "residue": [...], "granule_strength": ...}
for row in score_batch(
    [{"accession": "P1", "sequence": "MGGY..."}, {"accession": "P2", "sequence": "MASS..."}]
):
    print(row["accession"], row["single"])
```

- `single` — protein-level score. Default (distilled) route: the XGBoost
  surrogate's prediction of the PhaSePred web-archive catGRANULE value
  (absolute caliber). Formula routes (`route="paper"` / `"legacy-unaudited"`):
  the Z-normalized Equation-2 value. Pass `normalized=False` (or the CLI
  `--raw` flag) for the underlying raw Equation-2 / un-aligned values.
- `residue` — residue-level profile: an interior sliding-window mean of the
  per-residue Equation-1 scores. The default window is 51 residues
  (radius 25), covering positions `25..L-26` (1-indexed) — length `L-50`
  for `L >= 51`. This matches the PhaSePred web archive residue arrays.
  On the distilled route the profile is translated so its mean equals the
  surrogate `single` (per-protein affine alignment preserving the geometry).
- `granule_strength` — mean of the positive profile values (Equation 3).

## Command line

```console
$ catgranule predict --fasta proteins.fasta --out json
$ catgranule predict --fasta proteins.fasta --out csv --raw --route paper
$ catgranule check
```

Weight selection precedence: `--weights <path>` > `$CATGRANULE_WEIGHTS` >
`--route` (`distilled` default, `paper`, `legacy-unaudited`) >
packaged default.

## Weights and provenance

All scoring parameters live in replaceable JSON artifacts under
`src/catgranule/weights/` with a `provenance` block
(`route`: `distilled` | `paper` | `legacy-unaudited`, `source`, `date`,
optional `gate`):

- **distilled** (default, `mode="distilled"`): an XGBoost surrogate trained
  on the PhaSePred web-archive human eval-pool (2022-02-11 snapshot) that
  predicts the archive `single` directly from 181 per-sequence features
  (pure sequence + numpy, no external tools). Gate on the archive holdout:
  Spearman `rho`=0.998961, median `|delta|`=0.0177 — accepted with gap by
  user decision D21 (2026-09-07). The serialized model ships inside the
  package (`distilled_xgb_model.ubj`).
- **paper** (formula): the Bolognesi 2016 closed-form recipe — published
  coefficients + mmc1 Table S4 scales + yeast-proteome Z-normalization
  (recomputed n=5742, 2026-09-06). Alternate artifact (rho 0.9747 /
  delta 0.186, below gate); deterministic and tool-free.
- **legacy-unaudited** (formula): the pre-S4 reconstruction (same formula,
  older yeast normalization); kept for continuity and comparisons.

The redistribution terms of the upstream catGRANULE sources do not apply:
both artifacts are independent re-implementations of the published
equations / a distillation of the lab's own web archive.

## License

MIT.