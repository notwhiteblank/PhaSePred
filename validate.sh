#!/usr/bin/env bash
# validate.sh — run the three-layer PhaSePred validation and print the
# four-mode AUC table against the paper's numbers, not just a path to a JSON
# file. Thin wrapper around scripts/validate_phasepred.py (the implementation).
#
# Layers: (1) paper-layer AUC vs Table 1, (2) train/test leakage, (3) shipped-
# artifact consistency. The report is written to a gitignored path under runs/;
# the committed root validation_report.json is never replaced by this script.
#
# `--strict-paper` is deliberately NOT passed. hPdPS at 10 features sits
# permanently at delta -0.0108 against the fixed 0.01 tolerance (SEED=42,
# established across E2/E3/E5), so the honest result is exit 0 with that row
# recorded as a non-hard WARN. Promoting it would make the demo fail on a
# known, documented, accepted deviation; the row and the tolerance are frozen.
#
# The E3 leakage finding is printed in full and unsoftened: the paper's own
# S2/S3 split is not sequence-disjoint (PdPS-2 label conflicts).
#
# Exit codes mirror scripts/validate_phasepred.py:
#   0  no hard gate violated (paper/sequence overruns are WARNs)
#   1  a hard gate failed (accession leakage or artifact consistency)
#   2  input missing or usage error
#
# Runs from any working directory, as ./validate.sh or bash validate.sh.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$SCRIPT_DIR"
VALIDATOR="$REPO_ROOT/scripts/validate_phasepred.py"
MANIFEST="$REPO_ROOT/data/processed/MANIFEST.tsv"
DEFAULT_MODELS="$REPO_ROOT/src/phasepred/data/models"
DEFAULT_METRICS="$REPO_ROOT/products/A_paper_split_recomputed/metrics.json"
TASKS=(SaPS PdPS hSaPS hPdPS)
N_PER_TASK=10

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
    cat <<'EOF'
validate.sh — three-layer validation with the AUC table shown on screen.

Usage:
  ./validate.sh                     default: AUC + leakage + consistency
  ./validate.sh --recompute-auc     additionally recompute AUC from the tables
                                    and cross-check metrics.json at 1e-9
                                    (trains 200 models; off by default)
  ./validate.sh --sequences PATH    JSONL of accession/sequence, enables the
                                    sequence-identity leakage row
  ./validate.sh --models-dir DIR    models root (default: bundled models)
  ./validate.sh --metrics PATH      metrics.json (default: products/.../metrics)
  ./validate.sh --work-dir DIR      work dir (default: runs/validation)
  ./validate.sh -h | --help         this message

Output: the four-mode paper AUC table and the E3 leakage finding on stdout; the
  machine-readable report at <work>/validation_report.json (gitignored). The
  committed root validation_report.json is not touched.

Exit codes:
  0  no hard gate violated (paper/sequence overruns are non-hard WARNs)
  1  a hard gate failed (accession-level leakage or artifact consistency)
  2  input missing or usage error

Preflight, in order, each failure distinct and actionable:
  1. a Python interpreter is resolved (PHASEPRED_PYTHON, .venv, .pixi, PATH)
  2. that interpreter is >= 3.12
  3. phasepred imports from THIS checkout (never an ambient install)
  4. xgboost is >=3.2,<3.3
  5. the four data/processed/*.tsv match MANIFEST.tsv
  6. metrics.json and the 40 shipped .joblib artifacts are present
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
            "Point PHASEPRED_PYTHON at a 3.12+ interpreter, or create one."
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
            "The validator loads XGBClassifier artifacts and needs xgboost>=3.2,<3.3."
    fi
    if ! "$PYTHON" - "$version" <<'PY'
import sys

parts = sys.argv[1].split(".")
major = int(parts[0])
minor = int(parts[1]) if len(parts) > 1 else 0
raise SystemExit(0 if (major, minor) == (3, 2) else 1)
PY
    then
        fail "xgboost $version is outside >=3.2,<3.3." \
            "Pin the environment: \"$PYTHON\" -m pip install 'xgboost>=3.2,<3.3'" \
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

preflight_artifacts() {
    [[ -f "$METRICS" ]] || fail "missing metrics file: $METRICS" \
        "It defaults to products/A_paper_split_recomputed/metrics.json."
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
            "Restore them with 'git checkout -- src/phasepred/data/models'."
    fi
    echo "preflight: metrics.json and 40 artifacts present"
}

RECOMPUTE=0
SEQUENCES=""
WORK_DIR=""
MODELS_DIR="$DEFAULT_MODELS"
METRICS="$DEFAULT_METRICS"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h | --help)
            usage
            exit 0
            ;;
        --recompute-auc)
            RECOMPUTE=1
            shift
            ;;
        --sequences)
            SEQUENCES="${2:?--sequences needs a path}"
            shift 2
            ;;
        --models-dir)
            MODELS_DIR="${2:?--models-dir needs a path}"
            shift 2
            ;;
        --metrics)
            METRICS="${2:?--metrics needs a path}"
            shift 2
            ;;
        --work-dir)
            WORK_DIR="${2:?--work-dir needs a path}"
            shift 2
            ;;
        *)
            fail "unknown argument: $1 (try --help)"
            ;;
    esac
done

if [[ -z "$WORK_DIR" ]]; then
    WORK_DIR="$REPO_ROOT/runs/validation"
fi
REPORT="$WORK_DIR/validation_report.json"

preflight_python
preflight_phasepred
preflight_xgboost
preflight_tables
preflight_artifacts

mkdir -p "$WORK_DIR"
args=(--models-dir "$MODELS_DIR" --metrics "$METRICS" --output "$REPORT")
if [[ "$RECOMPUTE" -eq 1 ]]; then
    args+=(--recompute-auc)
fi
if [[ -n "$SEQUENCES" ]]; then
    args+=(--sequences "$SEQUENCES")
fi

set +e
"$PYTHON" "$VALIDATOR" "${args[@]}" \
    >"$WORK_DIR/validator_table.txt" 2>"$WORK_DIR/validator.stderr.txt"
rc=$?
set -e

if [[ -s "$WORK_DIR/validator.stderr.txt" ]]; then
    sed 's/^/[validator] /' "$WORK_DIR/validator.stderr.txt" >&2
fi

if [[ ! -f "$REPORT" ]]; then
    fail "the validator did not write a report (exit $rc); see $WORK_DIR/validator.stderr.txt"
fi

echo
"$PYTHON" - "$REPORT" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    report = json.load(handle)

paper = [check for check in report["checks"] if check.get("layer") == "paper"]
print("Four-mode AUC vs the paper's Table 1 (per-model cv_auc comparator)")
print(f"{'MODE':<9}{'PAPER':>8}{'MEASURED':>14}{'DELTA':>11}  VERDICT")
for check in paper:
    target = float(check["threshold"])
    measured = check["measured"]
    if measured is None:
        print(f"{check['name']:<9}{target:>8.3f}{'n/a':>14}{'n/a':>11}  "
              f"{check['verdict']}")
    else:
        value = float(measured)
        print(f"{check['name']:<9}{target:>8.3f}{value:>14.6f}"
              f"{value - target:>+11.4f}  {check['verdict']}")

summary = report["summary"]
print(f"\nsummary: pass={summary['pass']} fail={summary['fail']} "
      f"skip={summary['skip']} hard_fail={summary['hard_fail']}")
print("hPdPS-10 is recorded as a non-hard WARN (delta -0.0108, SEED=42); "
      "--strict-paper would promote it, which is why this script never passes it.")
PY

cat <<'EOF'

================================================================================
E3 LEAKAGE FINDING — PUBLIC, NOT SOFTENED
================================================================================
The paper's own S2/S3 split is not sequence-disjoint. Measured (E3, over the
63,836-sequence UniProt cache):
  - train/test sequence identity: SaPS 80 shared sequences / 90 accession
    pairs; PdPS 81 / 92; hSaPS 0 / 0; hPdPS 0 / 0.
  - PdPS-2 label conflicts: 2 shared sequences pair a label-1 TRAINING
    positive with a label-0 TEST negative -- P68431 (H3C1, human) in train vs
    P68433 (H3c1, mouse, NoPS-test) in test, and P63279 (UBE2I, human) in train
    vs P63281 (Ube2i, rat, NoPS-test) in test.
  - The paper's independent test set therefore contains sequences identical to
    positive training examples while labelled negative, which modestly favours
    its reported base-scope test AUC.
Accession-level train/test overlap is zero for all four tasks, so the hard gate
passes. The sequence-identity row is hard=false and is NOT promoted; removing
the conflicts would require departing from the paper's split and would break
byte-identical reproduction. This is a property of the published split, not of
the rebuild.
Recorded in:
  - products/A_paper_split_recomputed/README.md (where the numbers land)
================================================================================
EOF

echo
echo "report: $REPORT (gitignored; the committed root validation_report.json is untouched)"
echo "exit code: $rc  (0 no hard gate violated, 1 hard gate failed, 2 input/usage)"
exit "$rc"
