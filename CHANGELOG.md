# Changelog

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
