# Changelog

## [1.0.2] - 2026-09-28

Patch release for the 1.0.1 installation feedback. Model artifacts and
training data are unchanged.

### Fixed

- Use `python3` in the Ubuntu venv examples.
- Document DeepCoil's user-data installation path and show it in missing-tool
  guidance; keep the legacy checkout environment as a fallback.
- Show the PScore install command when `check-tools` reports it missing.
- Let the checkout PScore runner find DBS in the installed user-data copy and
  accept relative FASTA and output paths when invoked on its own.

## [1.0.1] - 2026-09-28

Patch release based on installation and prediction feedback. Model artifacts
and training data are unchanged.

### Fixed

- PScore's `install.sh --check` now checks the installed user-data copy instead
  of mistakenly checking for DBS files in the repository.
- PhaSePred sends DeepCoil smaller sequence groups to reduce peak memory use
  during batch prediction. A failed group reports a warning without discarding
  successful groups; sequences shorter than 20 residues get a specific warning.
- FCR remains available for sequences containing selenocysteine (`U`) even
  when the Hydropathy scale cannot score that residue.

### Improved

- `predict` identifies proteins with missing features before the model imputes
  those values and calculates a score.
- The READMEs clarify that `features-from-fasta` computes four features and
  that individual tool checks can be silent on success.

The local catGRANULE scale has no value for `U`, so catGRANULE remains missing
for such sequences and is imputed during prediction.
