#!/usr/bin/env bash
# train.sh — retrain all four PhaSePred models from the committed tables and
# compare the result with the 40 shipped model artifacts, byte for byte.
#
# This is a thin, auditable wrapper around scripts/train_phasepred.py (the
# implementation; docstrings and CLI flags live there). It never writes into
# src/phasepred/data/models/: it trains into a gitignored work directory and
# then compares against the shipped copy, so the byte-identity gate that seven
# phases depend on cannot be destroyed by running the demo.
#
# Reproduction is byte-exact only under xgboost >=3.2,<3.3. The preflight
# enforces that, so a wrong environment fails loudly instead of training 40
# non-identical files and producing a confusing diff.
#
# Exit codes:
#   0  every one of the 40 artifacts is byte-identical to the shipped copy
#   1  training succeeded but at least one artifact differs
#   2  preflight or usage failure (nothing was trained / nothing compared)
#
# Runs from any working directory, as ./train.sh or bash train.sh.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$SCRIPT_DIR"
TRAINER="$REPO_ROOT/scripts/train_phasepred.py"
MANIFEST="$REPO_ROOT/data/processed/MANIFEST.tsv"
DEFAULT_MODELS="$REPO_ROOT/src/phasepred/data/models"
TASKS=(SaPS PdPS hSaPS hPdPS)
N_PER_TASK=10

# Measured wall clock of a full default run on the reference environment
# (CPython 3.12, xgboost 3.2.0). Update if the reference timing changes.
FULL_RUN_WALL_CLOCK="~41 seconds (measured)"

fail() {
    printf 'error: %s\n' "$1" >&2
    shift
    local line
    for line in "$@"; do
        printf '  %s\n' "$line" >&2
    done
    exit 2
}

usage() {
    cat <<EOF
train.sh — one-command retraining from the committed data/processed/*.tsv.

Usage:
  ./train.sh                     full 4-task retrain, then byte-compare
  ./train.sh --smoke [N]         fast smoke retrain (N negatives/task, no CV)
  ./train.sh --compare-only      compare an existing work dir, do not train
  ./train.sh --work-dir DIR      work dir (default: runs/retrain[/-smoke])
  ./train.sh --models-dir DIR    shipped models to compare against
  ./train.sh -h | --help         this message

Modes:
  default        Trains SaPS/PdPS/hSaPS/hPdPS into <work>/models and writes
                 metrics into <work>/metrics, then compares every
                 8f_model_*.joblib against the 40 shipped artifacts.
  --smoke [N]    Wraps scripts/train_phasepred.py --smoke. A smoke run uses a
                 truncated negative pool, so its models are NOT byte-comparable;
                 it checks that 40 models were written and that the shipped
                 artifacts stayed unchanged instead.
  --compare-only Reads <work>/models and compares against --models-dir. This is
                 the non-vacuous re-check path (and how the comparison is
                 reverse-proved).

Output: runs/ (gitignored). The shipped models under
  src/phasepred/data/models/ are never written by any mode.

Exit codes:
  0  all 40 artifacts byte-identical to the shipped copy
  1  trained (or compared) but at least one artifact differs
  2  preflight or usage failure

Preflight, in order, each failure distinct and actionable:
  1. a Python interpreter is resolved (PHASEPRED_PYTHON, .venv, .pixi, PATH)
  2. that interpreter is >= 3.12
  3. phasepred imports from THIS checkout (never an ambient install)
  4. xgboost is >=3.2,<3.3 (3.4.x does not reproduce the published models)
  5. the four data/processed/*.tsv match MANIFEST.tsv
  6. the 40 shipped .joblib artifacts are present

Measured full-run wall clock on the reference environment: ${FULL_RUN_WALL_CLOCK}.
EOF
}

resolve_python() {
    local candidates=() candidate name found
    if [[ -n "${PHASEPRED_PYTHON:-}" ]]; then
        candidates+=("$PHASEPRED_PYTHON")
    fi
    candidates+=("$REPO_ROOT/.venv/bin/python" "$REPO_ROOT/.pixi/envs/default/bin/python")
    for name in python3.12 python3.13 python3 python; do
        found="$(command -v "$name" 2>/dev/null || true)"
        if [[ -n "$found" ]]; then
            candidates+=("$found")
        fi
    done
    for candidate in "${candidates[@]}"; do
        if [[ -x "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

preflight_python() {
    if ! PYTHON="$(resolve_python)"; then
        fail "no Python interpreter found." \
            "Set PHASEPRED_PYTHON, or provision this checkout with one of:" \
            "  pixi install" \
            "  python3.12 -m venv .venv && .venv/bin/pip install -e . -e packages/catgranule"
    fi
    local version
    version="$("$PYTHON" -c 'import sys; print(sys.version.split()[0])')"
    if ! "$PYTHON" -c \
        'import sys; raise SystemExit(0 if sys.version_info[:2] >= (3, 12) else 1)'; then
        fail "Python >= 3.12 required, but $PYTHON is $version." \
            "Point PHASEPRED_PYTHON at a 3.12+ interpreter, or create one:"
    fi
    echo "preflight: python $version ($PYTHON)"
}

preflight_phasepred() {
    local py_src imported
    py_src="$(cd "$REPO_ROOT/src" && pwd -P)"
    if ! imported="$(PYTHONPATH="$REPO_ROOT/src" "$PYTHON" -c \
        'import os, phasepred; print(os.path.realpath(phasepred.__file__))' 2>/dev/null)"; then
        fail "phasepred is not importable from $PYTHON." \
            "This checkout is not provisioned. Install it into the interpreter:" \
            "  \"$PYTHON\" -m pip install -e . -e packages/catgranule" \
            "(or use the pixi default environment)."
    fi
    case "$imported" in
        "$py_src"/*) ;;
        *)
            fail "phasepred resolved to $imported, not to this checkout ($py_src)." \
                "An ambient install would mask this tree's code. Provision a" \
                "dedicated environment and point PHASEPRED_PYTHON at it."
            ;;
    esac
    echo "preflight: phasepred from $imported"
}

preflight_xgboost() {
    local version
    if ! version="$("$PYTHON" -c 'import xgboost; print(xgboost.__version__)' 2>/dev/null)"; then
        fail "xgboost is not importable from $PYTHON." \
            "The trainer needs xgboost>=3.2,<3.3. Provision the environment:"
    fi
    if ! "$PYTHON" - "$version" <<'PY'
import sys

parts = sys.argv[1].split(".")
major = int(parts[0])
minor = int(parts[1]) if len(parts) > 1 else 0
raise SystemExit(0 if (major, minor) == (3, 2) else 1)
PY
    then
        fail "xgboost $version is outside the reproducible range >=3.2,<3.3." \
            "xgboost 3.4.x does NOT reproduce the published models: every one of" \
            "the 40 .joblib files would differ and the comparison would be" \
            "meaningless. Pin the environment instead of continuing:" \
            "  \"$PYTHON\" -m pip install 'xgboost>=3.2,<3.3'" \
            "(or use the pixi default environment, which pins ==3.2.0)."
    fi
    echo "preflight: xgboost $version"
}

preflight_tables() {
    [[ -f "$MANIFEST" ]] || fail "missing training-table manifest: $MANIFEST"
    local file scope split rows cols nbytes sha rest path actual
    while IFS=$'\t' read -r file scope split rows cols nbytes sha rest; do
        if [[ "$file" == "file" ]]; then
            continue
        fi
        path="$REPO_ROOT/data/processed/$file"
        [[ -f "$path" ]] || fail "missing training table: $path" \
            "The four committed tables under data/processed/ must be present."
        actual="$(sha256sum "$path" | awk '{print $1}')"
        if [[ "$actual" != "$sha" ]]; then
            fail "sha256 mismatch for data/processed/$file" \
                "  manifest: $sha" \
                "  on disk : $actual" \
                "The committed table is corrupted or has been edited. Restore it" \
                "with 'git checkout -- data/processed/$file'."
        fi
    done <"$MANIFEST"
    echo "preflight: 4 tables match MANIFEST.tsv"
}

preflight_models() {
    local task i missing=0
    for task in "${TASKS[@]}"; do
        for ((i = 0; i < N_PER_TASK; i++)); do
            if [[ ! -f "$MODELS_DIR/$task/8f_model_${i}.joblib" ]]; then
                missing=$((missing + 1))
            fi
        done
    done
    if ((missing > 0)); then
        fail "$missing of 40 shipped artifacts are missing under $MODELS_DIR." \
            "The shipped models are package data; restore them with" \
            "  git checkout -- src/phasepred/data/models"
    fi
    echo "preflight: 40 shipped artifacts present under $MODELS_DIR"
}

snapshot_shipped() {
    (
        cd "$DEFAULT_MODELS"
        find . -name '*.joblib' | sort | xargs sha256sum
    )
}

compare_models() {
    local identical=0 total=0 diffs=0 task i file out shipped a b
    local -a differing=()
    for task in "${TASKS[@]}"; do
        for ((i = 0; i < N_PER_TASK; i++)); do
            file="8f_model_${i}.joblib"
            out="$OUT_MODELS/$task/$file"
            shipped="$MODELS_DIR/$task/$file"
            total=$((total + 1))
            if [[ ! -f "$out" ]]; then
                differing+=("$task/$file (not produced at $out)")
                diffs=$((diffs + 1))
                continue
            fi
            a="$(sha256sum "$out" | awk '{print $1}')"
            b="$(sha256sum "$shipped" | awk '{print $1}')"
            if [[ "$a" == "$b" ]]; then
                identical=$((identical + 1))
                printf 'MATCH  %-28s %s\n' "$task/$file" "$a"
            else
                differing+=("$task/$file")
                diffs=$((diffs + 1))
                printf 'DIFF   %-28s\n' "$task/$file"
                printf '         work:    %s\n' "$a"
                printf '         shipped: %s\n' "$b"
            fi
        done
    done
    printf '\n%d/%d byte-identical\n' "$identical" "$total"
    if ((diffs > 0)); then
        printf 'differing artifacts:\n'
        local entry
        for entry in "${differing[@]}"; do
            printf '  - %s\n' "$entry"
        done
        return 1
    fi
    return 0
}

SMOKE_N=""
COMPARE_ONLY=0
WORK_DIR=""
MODELS_DIR="$DEFAULT_MODELS"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h | --help)
            usage
            exit 0
            ;;
        --smoke)
            if [[ "${2:-}" =~ ^[0-9]+$ ]]; then
                SMOKE_N="$2"
                shift 2
            else
                SMOKE_N="20"
                shift
            fi
            ;;
        --compare-only)
            COMPARE_ONLY=1
            shift
            ;;
        --work-dir)
            WORK_DIR="${2:?--work-dir needs a path}"
            shift 2
            ;;
        --models-dir)
            MODELS_DIR="${2:?--models-dir needs a path}"
            shift 2
            ;;
        *)
            fail "unknown argument: $1 (try --help)"
            ;;
    esac
done

if [[ -n "$SMOKE_N" && "$COMPARE_ONLY" -eq 1 ]]; then
    fail "--smoke and --compare-only are mutually exclusive"
fi

if [[ -z "$WORK_DIR" ]]; then
    if [[ -n "$SMOKE_N" ]]; then
        WORK_DIR="$REPO_ROOT/runs/retrain-smoke"
    else
        WORK_DIR="$REPO_ROOT/runs/retrain"
    fi
fi
OUT_MODELS="$WORK_DIR/models"
OUT_METRICS="$WORK_DIR/metrics"

# --- preflight -------------------------------------------------------------

preflight_python
preflight_phasepred
preflight_xgboost
preflight_tables
preflight_models

if [[ "$COMPARE_ONLY" -eq 1 ]]; then
    [[ -d "$OUT_MODELS" ]] || fail "no models to compare: $OUT_MODELS does not exist." \
        "Run ./train.sh first, or pass --work-dir pointing at a trained work dir."
    echo "compare-only: $OUT_MODELS vs $MODELS_DIR"
    compare_models
    exit $?
fi

before="$(snapshot_shipped)"

if [[ -n "$SMOKE_N" ]]; then
    echo "--- smoke retrain: train_phasepred.py --smoke $SMOKE_N ---"
    if ! "$PYTHON" "$TRAINER" --smoke "$SMOKE_N" \
        --output-dir "$OUT_MODELS" --metrics-dir "$OUT_METRICS"; then
        fail "smoke retrain failed under $PYTHON."
    fi
    written="$(find "$OUT_MODELS" -name '*.joblib' | wc -l)"
    ((written == 40)) || fail "smoke wrote $written artifacts, want 40" \
        "The smoke guard in scripts/train_phasepred.py refuses published paths;"
    for artifact in metrics.json manifest.json leakage_report.json; do
        [[ -f "$OUT_METRICS/$artifact" ]] \
            || fail "smoke did not write $artifact under $OUT_METRICS"
    done
    after="$(snapshot_shipped)"
    [[ "$after" == "$before" ]] \
        || fail "the shipped artifacts CHANGED during smoke retraining (should be impossible)"
    echo "smoke OK: 40 models under $OUT_MODELS; shipped artifacts unchanged"
    exit 0
fi

echo "--- full retrain: train_phasepred.py (all four tasks) ---"
if ! "$PYTHON" "$TRAINER" --output-dir "$OUT_MODELS" --metrics-dir "$OUT_METRICS"; then
    fail "full retraining failed under $PYTHON (see the traceback above)."
fi

written="$(find "$OUT_MODELS" -name '*.joblib' | wc -l)"
((written == 40)) || fail "training wrote $written artifacts, want 40" \
    "A partial run cannot be compared; inspect $OUT_MODELS."

after="$(snapshot_shipped)"
[[ "$after" == "$before" ]] \
    || fail "the shipped artifacts CHANGED during retraining (should be impossible)" \
        "src/phasepred/data/models/ must never be written by this script."

echo "--- comparison: $OUT_MODELS vs $MODELS_DIR ---"
if compare_models; then
    exit 0
fi
exit 1
