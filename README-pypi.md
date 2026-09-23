# PhaSePred

Predict the phase-separation propensity of a protein from its amino-acid sequence.
Implementation of [Chen et al. 2022 *PNAS*](https://doi.org/10.1073/pnas.2115369119).

## Install

```bash
pip install phasepred
```

Requires Python >= 3.12.

## Use

```bash
phasepred predict --fasta proteins.fasta --mode SaPS --output scores.csv
```

Four modes: `SaPS` and `PdPS` (any species, 8 features), `hSaPS` and `hPdPS`
(human, 10 features). The output CSV carries the score plus every feature value.

Input can also be UniProt accessions:

```bash
phasepred predict --ids "P35637,Q9Y2W1" --mode SaPS --output scores.csv
```

## Feature tools

Five of the nine feature components ship with the package and work immediately.
The rest need a one-time install:

```bash
phasepred check-tools
```

This prints a status row per component, and for every missing one the exact
command that installs it.

## More

Source, training data, retraining and AUC validation against the paper:
<https://github.com/notwhiteblank/PhaSePred>
