# PhaSePred

Predict the phase-separation propensity of a protein from its amino-acid sequence.
Implementation of [Chen et al. 2022 *PNAS*](https://doi.org/10.1073/pnas.2115369119).

## Install

```bash
pip install phasepred
```

Requires Python >= 3.12. The models and the DeepPhase table ship inside the
package.

Four feature components cannot be redistributed and need a one-time install from
the repository:

```bash
git clone https://github.com/notwhiteblank/PhaSePred.git && cd PhaSePred
bash tools/PScore/install.sh
bash tools/ESpritz/install.sh
bash tools/DeepCoil/install.sh
bash tools/PhosphoSitePlus/install.sh   # only needed for hSaPS / hPdPS
```

The installers put their data under `~/.local/share/phasepred/`, where the
package finds it automatically. `phasepred check-tools` prints a status row per
component, and for every missing one the exact command that installs it. A
complete install reports `9 OK`.

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

## More

Source, training data, retraining and AUC validation against the paper:
<https://github.com/notwhiteblank/PhaSePred>
