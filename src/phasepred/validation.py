"""Three-layer validation for the PhaSePred paper reproduction (E3).

The validator answers three independent questions about a checkout:

1. **AUC layer** -- does each shipped model's cross-validation AUC sit within
   0.01 of the paper's Table 1 target, measured with the *per-model* mean that
   the paper actually reports?  These rows are advisory (``hard=False``): the
   protocol fixes ``SEED=42`` and hPdPS-10 sits permanently at a delta of
   -0.0108, so a permanently red gate would be useless.  Callers that want the
   strict reading pass ``--strict-paper`` to :func:`exit_code_for`.
   :func:`recomputed_auc_layer` is the independent reading behind
   ``--recompute-auc``: it recomputes the fields from the committed tables and
   models and gates ``metrics.json`` at 1e-9.
2. **Leakage layer** -- is the train/test split disjoint per task?
   Accession overlap within one task is a hard failure: it is the check on our
   own pipeline, and it measures zero, as does the spec's label-1-across-the-split
   gate.  Identical sequences across the split are reported with ``hard=False``:
   the base-scope ``NoPS`` pool spans organisms, so cross-species orthologues can
   be sequence-identical (SaPS 80 sequences / 90 accession pairs, PdPS 81 / 92;
   the single-species human tasks show none), which is a property of the published
   Dataset S2/S3 splits.  A computed subset pairs label-1 training positives with
   label-0 test negatives (PdPS 2, SaPS 0), which modestly favours the paper's
   base-scope test AUC.  ``--strict-leakage`` promotes that row for callers who
   want it.  Cross-task accession sharing and shared test sheets are
   protocol-defined by the paper
   (``NoPS``/``hNoPS``/``PS-test``/``NoPS-test``/``hPS-test``/``hNoPS-test``
   each serve two tasks) and are reported without failing.  The sequence check
   skips -- never passes -- when sequences are unavailable, because the tables
   carry no sequences and the cache is gitignored.
3. **Artifact consistency layer** -- does each ``src/phasepred/data/models``
   manifest match the frozen protocol constants and the files on disk?

The report built by :func:`results_to_report` contains no wall-clock field:
inputs are identified by repository-relative path plus sha256, so two runs over
the same inputs produce byte-identical reports.

Every check carries a reverse proof in the tests: a plausible corruption of the
property it asserts makes it fail.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import xgboost
from sklearn.metrics import roc_auc_score

from phasepred.data import data_path, models_root, repo_root
from phasepred.models import (
    FEATURE_DEFINITIONS_VERSION,
    ModelManifestError,
    load_model_manifest,
    model_dir_is_stale,
)
from phasepred.training import (
    N_FOLDS,
    N_MODELS,
    NEG_RATIO,
    SEED,
    TASKS,
    ensemble_score,
    feature_columns,
    run_cross_validation,
    split_pos_neg,
)
from phasepred.training_data import (
    TASKS_BY_SCOPE,
    load_task_frame,
    sha256_file,
    table_path,
)

#: Paper Table 1 targets, `docs/RETRAIN_PROTOCOL.md` Sec. 5 (mineru_md:67).
#: The four shipped modes: base models at 8 features, human models at 10.
PAPER_TARGETS: dict[tuple[str, int], float] = {
    ("SaPS", 8): 0.862,
    ("PdPS", 8): 0.739,
    ("hSaPS", 10): 0.924,
    ("hPdPS", 10): 0.827,
}
#: Acceptance tolerance |delta AUC| <= 0.01. Never relax this.
PAPER_TOLERANCE = 0.01

#: Layer name for the independent AUC recomputation rows. Distinct from
#: ``"paper"`` so the recomputation is a hard artifact gate and the
#: ``--strict-paper`` switch does not also promote it.
AUC_RECOMPUTE_LAYER = "auc"
#: Recomputation must agree with ``metrics.json`` to 1e-9. Never relax this.
AUC_RECOMPUTE_TOLERANCE = 1e-9

#: Sentinel layer for the check row that carries input identity into the report.
INPUTS_LAYER = "inputs"

#: Name suffix identifying the sequence-identity report row within a task. The
#: row is reported but does not block (Ruling 2); ``--strict-leakage`` promotes
#: only rows matching this suffix.
SEQUENCE_ROW_SUFFIX = ": train/test sequence overlap"


class ValidationError(ValueError):
    """Raised when the validator cannot run (missing or unusable input)."""


@dataclass(frozen=True)
class CheckResult:
    """One row of validation output.

    ``verdict`` is ``pass``/``fail``/``skip``; ``hard`` marks the row as a
    gate.  Only ``hard and verdict == "fail"`` drives a non-zero exit code.
    """

    layer: str
    name: str
    verdict: str
    hard: bool
    measured: object
    threshold: object
    attribution: str


def auc_layer(
    metrics: Mapping[str, Mapping[str, Any]], *, strict_paper: bool = False
) -> list[CheckResult]:
    """Compare each paper row's per-model ``cv_auc`` against its Table 1 target.

    ``strict_paper`` is accepted for interface symmetry but deliberately does
    not alter the rows: the report always records the honest ``hard=False``,
    and the strict reading is applied by :func:`exit_code_for`.  The comparator
    is ``cv_auc`` (per-model mean), never ``cv_auc_ensemble`` -- the S5 errata.
    """
    del strict_paper
    results: list[CheckResult] = []
    for (task, n_features), target in PAPER_TARGETS.items():
        name = f"{task}-{n_features}"
        entry = metrics.get(task)
        measured = entry.get("cv_auc") if entry is not None else None
        if measured is None:
            results.append(
                CheckResult(
                    layer="paper",
                    name=name,
                    verdict="skip",
                    hard=False,
                    measured=None,
                    threshold=target,
                    attribution=(
                        f"metrics carries no cv_auc for task {task!r}; the AUC layer "
                        "cannot compare this row. metrics.json is produced by "
                        "scripts/train_phasepred.py."
                    ),
                )
            )
            continue
        value = float(measured)
        delta = value - target
        passed = abs(delta) <= PAPER_TOLERANCE
        attribution = ""
        if not passed:
            attribution = (
                "per-model cv_auc is the paper comparator (RETRAIN_PROTOCOL Sec. 5, "
                f"S5 errata); delta={delta:+.4f} exceeds the fixed tolerance. The "
                "protocol fixes SEED=42, so the miss is recorded rather than "
                "re-seeded or re-tuned; S5 already documents this boundary. The "
                "ensemble reading (cv_auc_ensemble) is not the paper comparator."
            )
        results.append(
            CheckResult(
                layer="paper",
                name=name,
                verdict="pass" if passed else "fail",
                hard=False,
                measured=value,
                threshold=target,
                attribution=attribution,
            )
        )
    return results


def _auc(labels: np.ndarray, scores: np.ndarray) -> float:
    """ROC AUC via ``sklearn.metrics.roc_auc_score``.

    Returns ``nan`` when the labels hold a single class, matching the trainer's
    finite-fold policy. The recomputation calls this directly rather than
    trusting the trainer's arithmetic; a unit test pins it against a
    hand-computed value.
    """
    if len(set(np.asarray(labels).tolist())) < 2:
        return float("nan")
    return float(roc_auc_score(labels, scores))


def _recompute_row(
    task: str, field: str, recomputed: float, recorded: float | None
) -> CheckResult:
    """One hard row comparing a recomputed AUC against its metrics field."""
    name = f"{task}: recomputed {field}"
    if recorded is None:
        return CheckResult(
            layer=AUC_RECOMPUTE_LAYER,
            name=name,
            verdict="fail",
            hard=True,
            measured={"recomputed": float(recomputed), "recorded": None},
            threshold=AUC_RECOMPUTE_TOLERANCE,
            attribution=(
                f"metrics.json carries no {field} for {task}; the committed metrics "
                "cannot be cross-checked against the committed tables and models. "
                "Regenerate metrics.json with scripts/train_phasepred.py."
            ),
        )
    delta = abs(float(recomputed) - float(recorded))
    passed = delta <= AUC_RECOMPUTE_TOLERANCE
    return CheckResult(
        layer=AUC_RECOMPUTE_LAYER,
        name=name,
        verdict="pass" if passed else "fail",
        hard=True,
        measured={
            "recomputed": float(recomputed),
            "recorded": float(recorded),
            "delta": delta,
        },
        threshold=AUC_RECOMPUTE_TOLERANCE,
        attribution=(
            ""
            if passed
            else (
                f"recomputed {field}={recomputed!r} differs from metrics.json "
                f"{recorded!r} by {delta:.3e}, above {AUC_RECOMPUTE_TOLERANCE:.0e}. The "
                "committed metrics no longer describe the committed tables and models; "
                "regenerate both with scripts/train_phasepred.py."
            )
        ),
    )


def recomputed_auc_layer(
    *,
    root: Path | None = None,
    models_dir: Path | None = None,
    metrics_path: Path | None = None,
) -> list[CheckResult]:
    """Recompute every AUC field from the committed tables and models.

    This is the independent reading of the AUC layer. It recomputes the
    per-model mean and ensemble CV AUC with
    :func:`phasepred.training.run_cross_validation` over the committed TSVs, and
    the S3 test AUC by loading ``models_dir/<task>/8f_model_*.joblib`` and
    scoring with :func:`phasepred.training.ensemble_score`. Each value is
    compared with the matching ``metrics.json`` field at
    ``AUC_RECOMPUTE_TOLERANCE`` and the row is a hard gate: drift means the
    committed metrics no longer describe the committed models.

    ``metrics_path`` defaults to the committed
    ``products/A_paper_split_recomputed/metrics.json``. ``run_all`` passes the
    file the caller selected with ``--metrics``, so a doctored metrics file is
    detected against the same input it replaces.

    The run trains 200 models and takes minutes, which is why the CLI gates it
    behind ``--recompute-auc`` (default off).

    :func:`~phasepred.training.run_cross_validation` prints one progress line
    per fold from inside its frozen protocol body, so the call is wrapped in
    ``contextlib.redirect_stdout(sys.stderr)`` (Ruling 8): this layer's consumer
    may be the ``--json`` report, and every diagnostic must stay off stdout.
    """
    metrics_file = (
        Path(metrics_path)
        if metrics_path is not None
        else data_path("products", "A_paper_split_recomputed", "metrics.json")
    )
    if not metrics_file.is_file():
        raise ValidationError(f"metrics file not found: {metrics_file}")
    metrics = json.loads(metrics_file.read_text(encoding="utf-8"))
    if not isinstance(metrics, dict):
        raise ValidationError(f"metrics file {metrics_file} is not a JSON object")
    model_root = Path(models_dir) if models_dir is not None else models_root()

    results: list[CheckResult] = []
    for task in TASKS:
        scope = "human" if task.startswith("h") else "base"
        train_frame = load_task_frame(scope, "train", task, root=root)
        test_frame = load_task_frame(scope, "test", task, root=root)
        cols = feature_columns(task)
        X_pos, X_neg_pool = split_pos_neg(train_frame, cols)
        with contextlib.redirect_stdout(sys.stderr):
            mean_auc, _fold_aucs, mean_ensemble, _fold_ensembles = run_cross_validation(
                X_pos, X_neg_pool, cols, SEED
            )

        directory = model_root / task
        files = sorted(directory.glob("8f_model_*.joblib"))
        if not files:
            raise ValidationError(
                f"no 8f_model_*.joblib under {directory}; regenerate the artifact with "
                "scripts/train_phasepred.py"
            )
        models = [joblib.load(path) for path in files]
        X_test = test_frame[cols].astype(float)
        y_test = test_frame["label"].astype(int).to_numpy()
        test_auc = _auc(y_test, ensemble_score(models, X_test))

        entry = metrics.get(task)
        entry = entry if isinstance(entry, dict) else {}
        recorded = {
            "cv_auc": entry.get("cv_auc"),
            "cv_auc_ensemble": entry.get("cv_auc_ensemble"),
            "test_auc_S3": entry.get("test_auc_S3"),
        }
        results.append(_recompute_row(task, "cv_auc", mean_auc, recorded["cv_auc"]))
        results.append(
            _recompute_row(
                task, "cv_auc_ensemble", mean_ensemble, recorded["cv_auc_ensemble"]
            )
        )
        results.append(
            _recompute_row(task, "test_auc_S3", test_auc, recorded["test_auc_S3"])
        )
    return results


def _accessions(frame: pd.DataFrame) -> set[str]:
    return {value.strip() for value in frame["UniprotEntry"].astype(str)}


def _accession_labels(frame: pd.DataFrame) -> dict[str, set[int]]:
    """Map each stripped accession to the set of integer labels on its rows."""
    labels: dict[str, set[int]] = {}
    for accession, label in zip(
        frame["UniprotEntry"].astype(str), frame["label"], strict=True
    ):
        try:
            value = int(label)
        except (TypeError, ValueError):
            continue
        labels.setdefault(accession.strip(), set()).add(value)
    return labels


def _accession_sheets(frame: pd.DataFrame) -> dict[str, set[str]]:
    """Map each stripped accession to the set of source sheets on its rows."""
    sheets: dict[str, set[str]] = {}
    for accession, sheet in zip(
        frame["UniprotEntry"].astype(str), frame["source_sheet"], strict=True
    ):
        sheets.setdefault(accession.strip(), set()).add(str(sheet))
    return sheets


def _accession_pairs(
    train_accessions: list[str],
    test_accessions: list[str],
    train_labels: Mapping[str, set[int]],
    test_labels: Mapping[str, set[int]],
    train_sheets: Mapping[str, set[str]],
    test_sheets: Mapping[str, set[str]],
) -> list[dict[str, object]]:
    """Both-sides accession pairs with their labels and source sheets."""
    pairs: list[dict[str, object]] = []
    for train_accession in train_accessions:
        for test_accession in test_accessions:
            pairs.append(
                {
                    "train_accession": train_accession,
                    "train_label": sorted(train_labels.get(train_accession, set())),
                    "train_sheet": sorted(train_sheets.get(train_accession, set())),
                    "test_accession": test_accession,
                    "test_label": sorted(test_labels.get(test_accession, set())),
                    "test_sheet": sorted(test_sheets.get(test_accession, set())),
                }
            )
    return pairs


def _shared_test_sheets(
    frames: Mapping[tuple[str, str], pd.DataFrame], tasks: Sequence[str]
) -> tuple[list[str], list[str]]:
    """Sheets shared by every task's test frame, and accessions on them.

    The paper's ``PS-test``/``NoPS-test``/``hPS-test``/``hNoPS-test`` sheets are
    reused across a scope by construction; identify them from ``source_sheet``
    rather than hard-coding accessions.
    """
    per_task: dict[str, set[str]] = {
        task: set(frames[(task, "test")]["source_sheet"].dropna().astype(str))
        for task in tasks
    }
    shared_sheets = sorted(set.intersection(*per_task.values())) if per_task else []
    if not shared_sheets:
        return [], []
    accessions: list[set[str]] = []
    for task in tasks:
        frame = frames[(task, "test")]
        mask = frame["source_sheet"].astype(str).isin(shared_sheets)
        accessions.append(set(frame.loc[mask, "UniprotEntry"].astype(str)))
    shared_accessions = sorted(set.intersection(*accessions)) if accessions else []
    return shared_sheets, shared_accessions


def leakage_layer(
    *, root: Path | None = None, sequences: Mapping[str, str] | None = None
) -> list[CheckResult]:
    """Accession- and sequence-level leakage checks over the training tables.

    For every ``(scope, task)``:

    * hard check (a): the train/test accession intersection within that task;
    * hard check (e): a label-1 accession that appears on both sides of the
      split -- the spec's positive-across-sheets gate, reported as its own row
      so the enumeration is visibly complete;
    * report check (b): the train/test sequence intersection within that task,
      computed only when ``sequences`` (accession -> sequence) is supplied and
      reported with ``hard=False`` (cross-species orthologues are inherent to
      the paper's base-scope split); otherwise the row is ``skip`` -- never
      ``pass``. The row carries both count readings (distinct shared sequences
      and accession pairs) and the label-conflict subset;
    * report check (c): accessions shared across the two tasks of a scope in
      their train frames (the ``NoPS``/``hNoPS`` protocol), ``hard=False``;
    * report check (d): test sheets shared across the two tasks of a scope,
      identified from ``source_sheet``, ``hard=False``.
    """
    results: list[CheckResult] = []
    for scope, tasks in TASKS_BY_SCOPE.items():
        frames: dict[tuple[str, str], pd.DataFrame] = {}
        for task in tasks:
            for split in ("train", "test"):
                frames[(task, split)] = load_task_frame(scope, split, task, root=root)

        for task in tasks:
            train_frame = frames[(task, "train")]
            test_frame = frames[(task, "test")]
            train_accessions = _accessions(train_frame)
            test_accessions = _accessions(test_frame)
            overlap = sorted(train_accessions & test_accessions)
            results.append(
                CheckResult(
                    layer="leakage",
                    name=f"{task}: train/test accession overlap",
                    verdict="fail" if overlap else "pass",
                    hard=True,
                    measured={"count": len(overlap), "overlap": overlap[:20]},
                    threshold=0,
                    attribution=(
                        f"{task} train and test share {len(overlap)} accession(s); the "
                        "paper protocol requires a disjoint split (RETRAIN_PROTOCOL "
                        "Sec. 8)."
                    )
                    if overlap
                    else "",
                )
            )

            train_labels = _accession_labels(train_frame)
            test_labels = _accession_labels(test_frame)
            train_sheets = _accession_sheets(train_frame)
            test_sheets = _accession_sheets(test_frame)
            positive_overlap = sorted(
                {a for a, labels in train_labels.items() if 1 in labels}
                & {a for a, labels in test_labels.items() if 1 in labels}
            )
            results.append(
                CheckResult(
                    layer="leakage",
                    name=f"{task}: positive in both train and test",
                    verdict="fail" if positive_overlap else "pass",
                    hard=True,
                    measured={
                        "count": len(positive_overlap),
                        "overlap": positive_overlap[:20],
                    },
                    threshold=0,
                    attribution=(
                        f"{task} carries {len(positive_overlap)} label-1 accession(s) on "
                        "both sides of the split (spec leakage check e); a positive in "
                        "both train and test inflates the test score and violates the "
                        "paper protocol (RETRAIN_PROTOCOL Sec. 8)."
                    )
                    if positive_overlap
                    else "",
                )
            )

            if sequences is None:
                results.append(
                    CheckResult(
                        layer="leakage",
                        name=f"{task}{SEQUENCE_ROW_SUFFIX}",
                        verdict="skip",
                        hard=False,
                        measured=None,
                        threshold=0,
                        attribution=(
                            "sequences unavailable: the committed tables carry no "
                            "sequences and data/interim/uniprot_cache.jsonl is "
                            "gitignored. Pass sequences=(accession -> sequence), or an "
                            "E5 fixture, to enable this report check."
                        ),
                    )
                )
            else:
                train_by_seq: dict[str, list[str]] = {}
                for accession in sorted(train_accessions):
                    sequence = sequences.get(accession)
                    if sequence is not None:
                        train_by_seq.setdefault(sequence, []).append(accession)
                test_by_seq: dict[str, list[str]] = {}
                for accession in sorted(test_accessions):
                    sequence = sequences.get(accession)
                    if sequence is not None:
                        test_by_seq.setdefault(sequence, []).append(accession)
                shared = sorted(set(train_by_seq) & set(test_by_seq))
                accession_pairs = 0
                examples: list[dict[str, object]] = []
                label_conflict_sequences: list[str] = []
                label_conflict_examples: list[dict[str, object]] = []
                for sequence in shared:
                    train_seq = train_by_seq[sequence]
                    test_seq = test_by_seq[sequence]
                    accession_pairs += len(train_seq) * len(test_seq)
                    if len(examples) < 5:
                        examples.extend(
                            _accession_pairs(
                                train_seq, test_seq, train_labels, test_labels,
                                train_sheets, test_sheets,
                            )
                        )
                    train_positive = any(
                        1 in train_labels.get(accession, set()) for accession in train_seq
                    )
                    test_negative = any(
                        0 in test_labels.get(accession, set()) for accession in test_seq
                    )
                    if train_positive and test_negative:
                        label_conflict_sequences.append(sequence)
                        if len(label_conflict_examples) < 5:
                            label_conflict_examples.extend(
                                _accession_pairs(
                                    train_seq, test_seq, train_labels, test_labels,
                                    train_sheets, test_sheets,
                                )
                            )
                conflict_clause = (
                    f"{len(label_conflict_sequences)} of those sequences pair a label-1 "
                    "training positive with a label-0 test negative: the paper's own "
                    "independent test set contains sequences identical to positive "
                    "training examples while labelled negative (a human positive against "
                    "a rodent orthologue labelled negative in test), which modestly "
                    "favours its reported base-scope test AUC. "
                    if label_conflict_sequences
                    else (
                        "None of those sequences pairs a label-1 training positive with "
                        "a label-0 test negative. "
                    )
                )
                results.append(
                    CheckResult(
                        layer="leakage",
                        name=f"{task}{SEQUENCE_ROW_SUFFIX}",
                        verdict="fail" if shared else "pass",
                        hard=False,
                        measured={
                            "shared_sequences": len(shared),
                            "accession_pairs": accession_pairs,
                            "label_conflict_sequences": len(label_conflict_sequences),
                            "examples": examples[:5],
                            "label_conflict_examples": label_conflict_examples[:5],
                        },
                        threshold=0,
                        attribution=(
                            f"{task} train and test share {len(shared)} identical "
                            f"sequence(s) across {accession_pairs} accession pair(s) under "
                            "different accessions. "
                            + conflict_clause
                            + "Cross-species orthologues account for the overlaps: the "
                            "base-scope NoPS negative pool spans organisms, so orthologues "
                            "can be identical, and the human tasks show none because they "
                            "are single-species. This is a property of the published "
                            "Dataset S2/S3 splits, not of this rebuild; removing the pairs "
                            "would require departing from the paper's split and would "
                            "break byte-identical reproduction. The row is reported "
                            "without blocking; accession-level overlap remains the hard "
                            "gate."
                        )
                        if shared
                        else "",
                    )
                )

        train_sets = [_accessions(frames[(task, "train")]) for task in tasks]
        shared_train = sorted(set.intersection(*train_sets)) if train_sets else []
        results.append(
            CheckResult(
                layer="leakage",
                name=f"{scope}: cross-task train accession sharing",
                verdict="pass",
                hard=False,
                measured={"count": len(shared_train), "accessions": shared_train[:20]},
                threshold=0,
                attribution=(
                    "the paper's NoPS/hNoPS negative sheets serve both tasks in a "
                    "scope, so shared accessions here are protocol-defined, not "
                    "leakage."
                ),
            )
        )

        shared_sheets, shared_test = _shared_test_sheets(frames, tasks)
        results.append(
            CheckResult(
                layer="leakage",
                name=f"{scope}: cross-task test accession sharing",
                verdict="pass",
                hard=False,
                measured={
                    "count": len(shared_test),
                    "accessions": shared_test[:20],
                    "sheets": shared_sheets,
                },
                threshold=0,
                attribution=(
                    "the paper's PS-test/NoPS-test/hPS-test/hNoPS-test sheets serve "
                    "both tasks in a scope (identified via source_sheet), so this "
                    "sharing is protocol-defined, not leakage."
                ),
            )
        )
    return results


def _verdict(
    name: str,
    ok: bool,
    *,
    measured: object,
    threshold: object,
    attribution: str,
    layer: str = "consistency",
    hard: bool = True,
) -> CheckResult:
    return CheckResult(
        layer=layer,
        name=name,
        verdict="pass" if ok else "fail",
        hard=hard,
        measured=measured,
        threshold=threshold,
        attribution="" if ok else attribution,
    )


def _major_minor(value: object) -> tuple[int, int] | None:
    parts = str(value).split(".")
    try:
        return (int(parts[0]), int(parts[1]))
    except (IndexError, ValueError):
        return None


def consistency_layer(models_dir: Path | None = None) -> list[CheckResult]:
    """Check each model directory's manifest against the frozen protocol.

    Every field is a hard gate except the xgboost version: a major/minor drift
    is an environment note, not an artifact defect, so it is reported with
    ``hard=False``.  A missing or unreadable manifest is a hard failure with an
    actionable remediation.  ``model_dir_is_stale`` is surfaced as its own row.
    """
    root = Path(models_dir) if models_dir is not None else models_root()
    runtime = xgboost.__version__
    results: list[CheckResult] = []
    for task in TASKS:
        directory = root / task
        try:
            manifest = load_model_manifest(directory)
        except ModelManifestError as exc:
            results.append(
                _verdict(
                    f"{task}: manifest",
                    False,
                    measured=str(directory),
                    threshold="readable manifest.json",
                    attribution=(
                        f"manifest is unreadable ({exc}); regenerate the artifact with "
                        "`python scripts/train_phasepred.py`."
                    ),
                )
            )
            continue

        if manifest is None:
            results.append(
                _verdict(
                    f"{task}: manifest",
                    False,
                    measured=None,
                    threshold="manifest.json",
                    attribution=(
                        f"{directory} has no manifest.json; regenerate the artifact "
                        "with `python scripts/train_phasepred.py` (RETRAIN_PROTOCOL "
                        "Sec. 6)."
                    ),
                )
            )
        else:
            raw = manifest.raw
            results.append(
                _verdict(
                    f"{task}: feature_definitions_version",
                    raw.get("feature_definitions_version") == FEATURE_DEFINITIONS_VERSION,
                    measured=raw.get("feature_definitions_version"),
                    threshold=FEATURE_DEFINITIONS_VERSION,
                    attribution=(
                        "manifest binds feature_definitions_version "
                        f"{raw.get('feature_definitions_version')!r} but the code "
                        f"implements {FEATURE_DEFINITIONS_VERSION!r}; retrain or "
                        "restore the v2022 artifact."
                    ),
                )
            )
            results.append(
                _verdict(
                    f"{task}: n_models",
                    raw.get("n_models") == N_MODELS,
                    measured=raw.get("n_models"),
                    threshold=N_MODELS,
                    attribution=(
                        f"manifest n_models={raw.get('n_models')!r}, expected {N_MODELS}."
                    ),
                )
            )
            files = sorted(directory.glob("8f_model_*.joblib"))
            declared = raw.get("n_models")
            expected_models: int = declared if isinstance(declared, int) else N_MODELS
            results.append(
                _verdict(
                    f"{task}: model_count",
                    len(files) == expected_models,
                    measured=len(files),
                    threshold=expected_models,
                    attribution=(
                        f"{directory} holds {len(files)} 8f_model_*.joblib file(s), "
                        f"manifest declares {expected_models}; regenerate with "
                        "`python scripts/train_phasepred.py`."
                    ),
                )
            )
            results.append(
                _verdict(
                    f"{task}: seed",
                    raw.get("seed") == SEED,
                    measured=raw.get("seed"),
                    threshold=SEED,
                    attribution=(
                        f"manifest seed={raw.get('seed')!r}, protocol fixes SEED={SEED}."
                    ),
                )
            )
            results.append(
                _verdict(
                    f"{task}: n_neg_ratio",
                    raw.get("n_neg_ratio") == NEG_RATIO,
                    measured=raw.get("n_neg_ratio"),
                    threshold=NEG_RATIO,
                    attribution=(
                        f"manifest n_neg_ratio={raw.get('n_neg_ratio')!r}, expected "
                        f"{NEG_RATIO}."
                    ),
                )
            )
            results.append(
                _verdict(
                    f"{task}: n_folds",
                    raw.get("n_folds") == N_FOLDS,
                    measured=raw.get("n_folds"),
                    threshold=N_FOLDS,
                    attribution=(
                        f"manifest n_folds={raw.get('n_folds')!r}, expected {N_FOLDS}."
                    ),
                )
            )
            expected_columns = feature_columns(task)
            results.append(
                _verdict(
                    f"{task}: feature_columns",
                    list(raw.get("feature_columns") or []) == expected_columns,
                    measured=list(raw.get("feature_columns") or []),
                    threshold=expected_columns,
                    attribution=(
                        "manifest feature_columns differ from "
                        f"phasepred.training.feature_columns({task!r}); the artifact "
                        "binds a different feature set."
                    ),
                )
            )
            results.append(
                _verdict(
                    f"{task}: xgboost_version",
                    _major_minor(raw.get("xgboost_version")) == _major_minor(runtime),
                    measured=raw.get("xgboost_version"),
                    threshold=runtime,
                    attribution=(
                        f"manifest xgboost_version={raw.get('xgboost_version')!r}, "
                        f"runtime is {runtime!r}; a major/minor difference is an "
                        "environment note, not an artifact defect."
                    ),
                    hard=False,
                )
            )

        stale, reason = model_dir_is_stale(directory)
        results.append(
            _verdict(
                f"{task}: staleness",
                not stale,
                measured=stale,
                threshold=False,
                attribution=reason or "model directory binds the current feature definitions.",
            )
        )
    return results


def _relative(path: Path) -> str:
    """Repository-relative posix path, or a non-path marker when outside it.

    The fallback deliberately never returns an absolute path: committed
    reports and manifests must not leak the machine's directory layout
    (``/mnt/data``, ``/home/zrc``).  Outside inputs are identified by file name
    plus an explicit marker.
    """
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(repo_root()).as_posix()
    except ValueError:
        return f"{resolved.name} (outside repository)"


def _mapping_sha256(mapping: Mapping[str, str]) -> str:
    digest = hashlib.sha256()
    for key in sorted(mapping):
        digest.update(str(key).encode("utf-8"))
        digest.update(b"\0")
        digest.update(str(mapping[key]).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def _inputs_result(
    *, root: Path | None, metrics_path: Path, sequences: Mapping[str, str] | None
) -> CheckResult:
    tables: dict[str, object] = {}
    for scope in TASKS_BY_SCOPE:
        for split in ("train", "test"):
            path = table_path(scope, split, root=root)
            tables[f"{scope}/{split}"] = {
                "path": _relative(path),
                "sha256": sha256_file(path) if path.is_file() else None,
            }
    payload: dict[str, object] = {
        "tables": tables,
        "metrics": {"path": _relative(metrics_path), "sha256": sha256_file(metrics_path)},
    }
    if sequences is not None:
        payload["sequences"] = {
            "count": len(sequences),
            "sha256": _mapping_sha256(sequences),
        }
    return CheckResult(
        layer=INPUTS_LAYER,
        name=INPUTS_LAYER,
        verdict="pass",
        hard=False,
        measured=payload,
        threshold=None,
        attribution="",
    )


def run_all(
    *,
    root: Path | None = None,
    models_dir: Path | None = None,
    metrics_path: Path | None = None,
    sequences: Mapping[str, str] | None = None,
    strict_paper: bool = False,
    recompute_auc: bool = False,
) -> list[CheckResult]:
    """Run all three layers and append the input-identity carrier row.

    ``metrics_path`` defaults to the committed
    ``products/A_paper_split_recomputed/metrics.json``.  A missing or non-object
    metrics file raises :class:`ValidationError`; the CLI maps that to exit 2.

    ``recompute_auc`` is off by default.  When set, the AUC fields are
    independently recomputed from the committed tables and models and compared
    with ``metrics_path`` at 1e-9 (hard gate).  The recomputation trains 200
    models, so the CLI exposes it as ``--recompute-auc``.
    """
    metrics_file = (
        Path(metrics_path)
        if metrics_path is not None
        else data_path("products", "A_paper_split_recomputed", "metrics.json")
    )
    if not metrics_file.is_file():
        raise ValidationError(f"metrics file not found: {metrics_file}")
    metrics = json.loads(metrics_file.read_text(encoding="utf-8"))
    if not isinstance(metrics, dict):
        raise ValidationError(f"metrics file {metrics_file} is not a JSON object")

    results: list[CheckResult] = []
    results.extend(auc_layer(metrics, strict_paper=strict_paper))
    results.extend(leakage_layer(root=root, sequences=sequences))
    results.extend(consistency_layer(models_dir=models_dir))
    if recompute_auc:
        results.extend(
            recomputed_auc_layer(
                root=root, models_dir=models_dir, metrics_path=metrics_file
            )
        )
    results.append(_inputs_result(root=root, metrics_path=metrics_file, sequences=sequences))
    return results


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        return [_jsonable(item) for item in sorted(value, key=str)]
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def results_to_report(results: Sequence[CheckResult]) -> dict[str, object]:
    """Serialize results into the stable, clock-free report schema.

    The input-identity carrier row (``layer == "inputs"``) becomes the top-level
    ``inputs`` map and is not repeated in ``checks``.  ``generated_from`` names
    the producing module.  No key in the output is a wall-clock reading.
    """
    checks: list[dict[str, object]] = []
    inputs: dict[str, object] = {}
    for result in results:
        if result.layer == INPUTS_LAYER:
            if isinstance(result.measured, Mapping):
                inputs = dict(result.measured)
            continue
        checks.append(
            {
                "layer": result.layer,
                "name": result.name,
                "verdict": result.verdict,
                "hard": bool(result.hard),
                "measured": _jsonable(result.measured),
                "threshold": _jsonable(result.threshold),
                "attribution": result.attribution,
            }
        )
    summary = {
        "pass": sum(1 for check in checks if check["verdict"] == "pass"),
        "fail": sum(1 for check in checks if check["verdict"] == "fail"),
        "skip": sum(1 for check in checks if check["verdict"] == "skip"),
        "hard_fail": sum(
            1 for check in checks if check["verdict"] == "fail" and check["hard"]
        ),
    }
    return {
        "generated_from": "phasepred.validation",
        "checks": checks,
        "summary": summary,
        "thresholds": {
            "paper_tolerance": PAPER_TOLERANCE,
            "paper_targets": {
                f"{task}-{n_features}": target
                for (task, n_features), target in PAPER_TARGETS.items()
            },
        },
        "inputs": inputs,
    }


def exit_code_for(
    results: Sequence[CheckResult],
    *,
    strict_paper: bool = False,
    strict_leakage: bool = False,
) -> int:
    """Tiered exit code (Ruling 1).

    ``1`` on any ``hard and verdict == "fail"``.  Otherwise the two strict
    switches each promote only their own layer: ``strict_paper`` promotes a
    ``layer == "paper"`` failure, ``strict_leakage`` promotes the sequence-
    identity report row.  Neither promotes the consistency layer's xgboost
    version-drift note, which is an environment observation.  ``skip`` never
    drives the exit code.
    """
    if any(result.hard and result.verdict == "fail" for result in results):
        return 1
    if strict_paper and any(
        result.layer == "paper" and result.verdict == "fail" for result in results
    ):
        return 1
    if strict_leakage and any(
        result.verdict == "fail"
        and result.layer == "leakage"
        and result.name.endswith(SEQUENCE_ROW_SUFFIX)
        for result in results
    ):
        return 1
    return 0


def _cell(value: object) -> str:
    if value is None:
        return "-"
    if isinstance(value, str):
        return value
    return json.dumps(_jsonable(value), sort_keys=True)


def format_report_table(results: Sequence[CheckResult]) -> str:
    """Render a fixed-width text table; FAIL rows carry an indented attribution."""
    header = (
        f"{'LAYER':<10} {'NAME':<40} {'VERDICT':<7} {'HARD':<5} {'MEASURED':<24} THRESHOLD"
    )
    lines = [header, "-" * len(header)]
    for result in results:
        if result.layer == INPUTS_LAYER:
            continue
        lines.append(
            f"{result.layer:<10} {result.name:<40} {result.verdict:<7} "
            f"{str(result.hard):<5} {_cell(result.measured):<24} {_cell(result.threshold)}"
        )
        if result.verdict == "fail" and result.attribution:
            lines.append(f"    {result.attribution}")
    return "\n".join(lines)
