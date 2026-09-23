"""Manifest-based external tool registry and path resolution.

Since S3 every feature tool is a self-contained package under
``tools/<TOOL>/`` with a machine-readable ``manifest.toml`` (schema in
``docs/FEATURES.md`` 附注 A and ``docs/TOOL_LICENSES.md`` §12). This module
scans those manifests to build a registry and resolves each tool with the
cascade::

    PHASEPRED_<TOOL>_DIR  →  user data dir / checkout  →  tools/<TOOL>/  →  PATH  →  error

The user-data tier (added in E4) probes ``user_data_root()/<tool>/<entry>``
through :func:`phasepred.data.resolve_optional_data`, so a complete contract
package unpacked under ``$PHASEPRED_DATA_ROOT`` / ``$XDG_DATA_HOME/phasepred``
/ ``~/.local/share/phasepred`` is preferred over the checkout copy. Components
without an entry in ``data._COMPONENT_LEGACY_NAMES`` have no user tier. Per
Ruling 9 the tier is reported honestly: user-controlled roots as ``user`` and
the checkout fallback as ``repo``.

``check_all_tools()`` returns a ``ToolStatus`` list (carrying the manifest's
``license`` / ``redistributable`` for the ``phasepred check-tools`` table),
verifies *runtime usability* through each package's ``./run --check``, and
is used by predict()'s preflight gate. catGRANULE is resolved as a package
import (``source=package``), not a tool contract. DeepPhase and
PhosphoSitePlus are data rows appended synthetically (E4: the registry is 9
rows covering the 10 feature columns; LocalCIDER supplies two columns).

The resolver never invokes a tool for its core run — package ``run``
scripts are the single execution entrypoint (see ``phasepred.tools``).

``ToolStatus.source`` comes from :data:`SOURCE_VOCABULARY`; both
``format_status_table`` and ``statuses_to_json`` canonicalise through the shared
:func:`_render_source` table, so they render exactly the same vocabulary
(``env`` is spelled ``env:<VAR>``).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tomllib
from dataclasses import dataclass, replace
from pathlib import Path

from phasepred.data import (
    bundled_path,
    data_path,
    repo_root,
    resolve_optional_data,
    user_data_root,
)

TOOLS_ROOT: Path = data_path("tools")


def tools_root() -> Path:
    """Resolve the tool-package scan root (three-stage fallback, Ruling 6).

    Preference, matching :func:`phasepred.data.models_root`'s style:

    1. ``$PHASEPRED_DATA_ROOT/tools`` when the override is set and that
       directory exists (external data root).
    2. the wheel-bundled contracts (``src/phasepred/data/tools``) when the
       override has no ``tools/`` subtree. Ruling 7: this tier was measured
       safe to scan — six validated manifests, no ``ManifestSchemaError`` —
       so a wheel user who follows ``docs/INSTALL.md`` and sets
       ``PHASEPRED_DATA_ROOT`` keeps the full 9-row registry instead of
       collapsing to the three synthetic rows (the pre-E5 D42 defect).
    3. the dev checkout (``repo_root()/tools``, full entities) as the last
       resort.

    Without ``PHASEPRED_DATA_ROOT`` the checkout ``tools/`` is preferred over
    the bundled copy (in a dev checkout the repository contracts are the
    canonical entities; the bundled copies are a lean wheel fallback). Never
    creates directories, and evaluates the environment on every call.
    """
    override = os.environ.get("PHASEPRED_DATA_ROOT")
    if override:
        override_tools = Path(override).expanduser().resolve() / "tools"
        if override_tools.is_dir():
            return override_tools
        bundled = bundled_path("data", "tools")
        if bundled.is_dir():
            return bundled
        return repo_root() / "tools"
    repo = data_path("tools")
    if repo.is_dir():
        return repo
    return bundled_path("data", "tools")

KIND_ENUM = {
    "binary",
    "java",
    "perl",
    "python",
    "env-prefix",
    "python-library",
    # Ruling 11: data-only contracts (E4: tools/PhosphoSitePlus) pin a download
    # rather than a runnable tool. ``load_manifests`` validates them but excludes
    # them from the tool registry, so their feature column keeps its single
    # synthetic row and the registry stays 9 rows.
    "data",
}

DISPLAY_NAMES = {
    "seg": "SEG",
    "plaac": "PLAAC",
    "pscore": "PScore",
    "espritz": "ESpritz",
    "deepcoil": "DeepCoil",
    "localcider": "LocalCIDER",
}

# Preferred table order (matches the S1 tool-contract draft order).
_ORDER = {key: i for i, key in enumerate(DISPLAY_NAMES)}

# Feature columns each contract accounts for (FEATURES.md). The 9 registry rows
# together cover all 10 feature columns; LocalCIDER supplies two native columns.
FEATURE_COLUMNS: dict[str, tuple[str, ...]] = {
    "seg": ("LCR",),
    "plaac": ("PLAAC",),
    "pscore": ("PScore",),
    "espritz": ("IDR",),
    "deepcoil": ("DeepCoil",),
    "localcider": ("Hydropathy", "FCR"),
}

# Registered ``ToolStatus.source`` values. Both ``format_status_table`` and
# ``statuses_to_json`` canonicalise through :func:`_render_source`, so they can
# only emit these values; "env" is written "env:<VAR>".
SOURCE_VOCABULARY = ("env", "user", "repo", "vendored", "path", "package", "missing")

# Ruling 9: one tier → source table. ``data.resolve_optional_data`` reports the
# tier that supplied a file (``override``/``xdg``/``home``/``repo``/``none``);
# this is the only place those become the public ``source`` value. The non-tier
# rows are identity entries so that the same single table also backs
# :func:`_render_source`, which both renderers consult — a resolver and a
# renderer therefore cannot disagree about provenance.
_TIER_TO_SOURCE: dict[str, str] = {
    "override": "user",
    "xdg": "user",
    "home": "user",
    "repo": "repo",
    "none": "user",
    "env": "env",
    "user": "user",
    "vendored": "vendored",
    "path": "path",
    "package": "package",
    "missing": "missing",
}


def _source_for_tier(tier: str) -> str:
    """Map a :func:`phasepred.data.resolve_optional_data` tier to a source value."""
    try:
        return _TIER_TO_SOURCE[tier]
    except KeyError as e:
        raise ValueError(f"unregistered optional-data tier {tier!r}") from e


def _render_source(source: str) -> str:
    """Canonicalise a source through the shared mapping table (Ruling 9).

    Both :func:`format_status_table` and :func:`statuses_to_json` call this, so
    the table and :data:`SOURCE_VOCABULARY` are load-bearing rather than a test
    oracle. ``env:<VAR>`` keeps its variable name; any other value must be
    registered or the renderer fails loudly.
    """
    base = source.split(":", 1)[0]
    if base not in _TIER_TO_SOURCE or base not in SOURCE_VOCABULARY:
        raise ValueError(f"unregistered tool source {source!r}")
    return source if ":" in source else _TIER_TO_SOURCE[base]


class PhaSePredToolNotFound(RuntimeError):
    """Raised when a required tool cannot be located by the resolver."""


class ManifestSchemaError(ValueError):
    """Raised when a ``tools/<TOOL>/manifest.toml`` violates the schema."""


@dataclass(frozen=True)
class ToolManifest:
    """Parsed and validated contents of a tool-package ``manifest.toml``."""

    key: str
    kind: str
    entry: str
    license: str
    redistributable: bool
    runtime: dict[str, str]
    path_executables: tuple[str, ...]
    data: tuple[dict[str, str], ...]
    output_format: str
    output_description: str
    install_command: str
    install_hint: str
    package_dir: Path

    @property
    def env_var(self) -> str:
        return f"PHASEPRED_{self.key.upper()}_DIR"

    @property
    def entry_path(self) -> Path:
        return self.package_dir / self.entry


@dataclass
class ToolStatus:
    """Resolution outcome for a single tool."""

    name: str
    ok: bool
    source: str  # "env:<VAR>", "user", "vendored", "path", "package", or "missing"
    path: Path | None
    required_for_predict: bool
    license: str = ""
    redistributable: bool = True
    kind: str = ""
    install_hint: str = ""
    feature_columns: tuple[str, ...] = ()


# ---------------------------------------------------------------------------
# Manifest loading + schema validation
# ---------------------------------------------------------------------------


def _display(key: str) -> str:
    return DISPLAY_NAMES.get(key, key.upper())


def _require(v: object, field: str, where: str) -> object:
    if v is None or v == "":
        raise ManifestSchemaError(f"manifest {where}: missing required field '{field}'")
    return v


def _parse_manifest(package_dir: Path) -> ToolManifest:
    toml_path = package_dir / "manifest.toml"
    root = tools_root()
    try:
        if toml_path.is_relative_to(root):
            where = str(toml_path.relative_to(root))
        else:
            where = str(toml_path)
    except ValueError:
        where = str(toml_path)
    try:
        raw = tomllib.loads(toml_path.read_text(encoding="utf-8"))
    except tomllib.TOMLDecodeError as e:
        raise ManifestSchemaError(f"manifest {where}: invalid TOML: {e}") from e

    key = str(_require(raw.get("name"), "name", where)).lower()
    kind = str(_require(raw.get("kind"), "kind", where))
    if kind not in KIND_ENUM:
        raise ManifestSchemaError(
            f"manifest {where}: kind '{kind}' not in {sorted(KIND_ENUM)}"
        )
    entry = str(_require(raw.get("entry"), "entry", where))
    license_ = str(_require(raw.get("license"), "license", where))
    redistributable = raw.get("redistributable")
    if not isinstance(redistributable, bool):
        raise ManifestSchemaError(
            f"manifest {where}: 'redistributable' must be a boolean"
        )

    runtime_raw = raw.get("runtime", {})
    if not isinstance(runtime_raw, dict):
        raise ManifestSchemaError(f"manifest {where}: 'runtime' must be a table")
    runtime = {str(k): str(v) for k, v in runtime_raw.items()}

    path_exes_raw = raw.get("path_executables", [])
    if not isinstance(path_exes_raw, list) or not all(
        isinstance(p, str) for p in path_exes_raw
    ):
        raise ManifestSchemaError(f"manifest {where}: 'path_executables' must be a list of strings")
    path_executables = tuple(path_exes_raw)

    data_raw = raw.get("data", [])
    data: list[dict[str, str]] = []
    if data_raw:
        if not isinstance(data_raw, list):
            raise ManifestSchemaError(f"manifest {where}: 'data' must be an array of tables")
        for item in data_raw:
            if not isinstance(item, dict) or "path" not in item:
                raise ManifestSchemaError(
                    f"manifest {where}: each [[data]] entry needs a 'path' string"
                )
            data.append(
                {
                    "path": str(item["path"]),
                    "note": str(item.get("note", "")),
                }
            )

    output = raw.get("output", {})
    if not isinstance(output, dict):
        raise ManifestSchemaError(f"manifest {where}: 'output' must be a table")
    output_format = str(_require(output.get("format"), "output.format", where))
    output_description = str(output.get("description", ""))

    install = raw.get("install", {})
    if not isinstance(install, dict):
        raise ManifestSchemaError(f"manifest {where}: 'install' must be a table")
    install_command = str(_require(install.get("command"), "install.command", where))
    install_hint = str(install.get("hint", ""))

    return ToolManifest(
        key=key,
        kind=kind,
        entry=entry,
        license=license_,
        redistributable=redistributable,
        runtime=runtime,
        path_executables=path_executables,
        data=tuple(data),
        output_format=output_format,
        output_description=output_description,
        install_command=install_command,
        install_hint=install_hint,
        package_dir=package_dir,
    )


def load_manifests(tools_root_path: Path | None = None) -> list[ToolManifest]:
    """Scan ``tools/*/manifest.toml`` and return validated manifests, sorted.

    A missing tools root yields ``[]`` (no crash) so a pip-installed wheel
    without a checkout degrades gracefully; a parse/validation failure raises
    ``ManifestSchemaError`` so the failure is loud rather than silently
    skipping a tool.

    ``kind = "data"`` contracts (Ruling 11) are parsed and validated, then
    excluded: they are not runnable registry rows. E4's PhosphoSitePlus contract
    is reported through the synthetic :func:`_phosphosite_status` row instead,
    which keeps the registry at 9 rows and its ``required_for_predict`` false.
    """
    root = tools_root_path or tools_root()
    if not root.is_dir():
        return []
    manifests: list[ToolManifest] = []
    for package_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        if (package_dir / "manifest.toml").is_file():
            manifest = _parse_manifest(package_dir)
            if manifest.kind == "data":
                continue
            manifests.append(manifest)
    return sorted(manifests, key=lambda m: _ORDER.get(m.key, len(_ORDER)))


def _find_manifest(key: str) -> ToolManifest:
    for manifest in load_manifests():
        if manifest.key == key:
            return manifest
    raise KeyError(f"no tool manifest registered for '{key}'")


# ---------------------------------------------------------------------------
# Resolution (env -> vendored package dir -> PATH -> missing)
# ---------------------------------------------------------------------------


def _status_for(
    manifest: ToolManifest,
    *,
    ok: bool,
    source: str,
    path: Path | None,
) -> ToolStatus:
    return ToolStatus(
        name=_display(manifest.key),
        ok=ok,
        source=source,
        path=path,
        required_for_predict=True,
        license=manifest.license,
        redistributable=manifest.redistributable,
        kind=manifest.kind,
        install_hint=manifest.install_hint,
        feature_columns=FEATURE_COLUMNS.get(manifest.key, ()),
    )


def _env_prefix_ready(manifest: ToolManifest, prefix: Path) -> bool:
    """Whether ``prefix`` holds a usable isolated env for an ``env-prefix`` tool.

    Mirrors ``bin/<key>`` / ``bin/python`` executability in
    ``tools/<TOOL>/run`` and ``tools/<TOOL>/install.sh``. Those two shell
    scripts and this function encode the same Ruling 10 tier order; a test
    (``test_deepcoil_three_way_prefix_agreement``) fails if they drift.
    """
    return os.access(prefix / "bin" / manifest.key, os.X_OK) or os.access(
        prefix / "bin" / "python", os.X_OK
    )


def _resolve_env_prefix(manifest: ToolManifest) -> ToolStatus | None:
    """Resolve an ``env-prefix`` contract (DeepCoil) to its vendored runner.

    The entity is an isolated conda env, not a file ``entry``, so the generic
    entry probe cannot see it. Tier order mirrors ``tools/DeepCoil/run`` and
    ``tools/DeepCoil/install.sh`` (Ruling 10):

    1. ``$<KEY>_ENV_PREFIX`` — authoritative when set, even if it points nowhere.
    2. ``user_data_root()/envs/<key>`` — the installer's write target
       (``source="user"``).
    3. ``repo_root()/<runtime.env_prefix>`` — the legacy checkout env
       (``source="vendored"``: a local build artifact, not user data).

    The resolved *runner* is always the vendored contract ``run``; ``source``
    reports which tier supplied the env. Returns ``None`` when no tier holds a
    usable prefix, and the caller reports MISSING.
    """
    env_var = f"{manifest.key.upper()}_ENV_PREFIX"
    env_val = os.environ.get(env_var)
    if env_val:
        if _env_prefix_ready(manifest, Path(env_val).expanduser()):
            return _status_for(
                manifest, ok=True, source=f"env:{env_var}", path=manifest.entry_path
            )
        return None

    user_prefix = user_data_root() / "envs" / manifest.key
    if _env_prefix_ready(manifest, user_prefix):
        return _status_for(manifest, ok=True, source="user", path=manifest.entry_path)

    legacy_rel = manifest.runtime.get("env_prefix")
    if legacy_rel and _env_prefix_ready(manifest, repo_root() / legacy_rel):
        return _status_for(
            manifest, ok=True, source="vendored", path=manifest.entry_path
        )
    return None


def _resolve(manifest: ToolManifest) -> ToolStatus:
    """Run the 5-tier cascade for one tool (without usability probe).

    Order: ``PHASEPRED_<TOOL>_DIR`` → user data dir / checkout → vendored
    ``tools/<TOOL>/`` → ``PATH`` (manifest ``path_executables``) → missing.
    The user-data tier goes through
    :func:`phasepred.data.resolve_optional_data` and only exists for components
    registered there; a tier-3 checkout hit reports ``source="repo"``
    (Ruling 9), not ``"user"``.

    ``kind = "env-prefix"`` contracts (DeepCoil) are the exception: their entity
    is an isolated env rather than a file, so they resolve through
    :func:`_resolve_env_prefix`, whose tier order matches the contract's ``run``
    and ``install.sh``.
    """

    # 1. PHASEPRED_<TOOL>_DIR
    env_val = os.environ.get(manifest.env_var)
    if env_val:
        cand = Path(env_val).expanduser()
        entry = cand / manifest.entry
        if entry.is_file():
            return _status_for(manifest, ok=True, source=f"env:{manifest.env_var}", path=entry)

    # 1b. env-prefix contracts: resolve the isolated env, not the entry file.
    if manifest.kind == "env-prefix":
        resolved = _resolve_env_prefix(manifest)
        if resolved is not None:
            return resolved
        return _status_for(manifest, ok=False, source="missing", path=None)

    # 2. user data dir / checkout fallback: a complete package unpacked under
    #    the user-writable root, or a tier-3 checkout hit. resolve_optional_data
    #    raises KeyError for components with no registered user layout, which
    #    simply means "no user tier". Ruling 9: report the resolved tier
    #    honestly ("repo" for a checkout hit) instead of labelling every hit
    #    "user".
    try:
        user_entry, tier = resolve_optional_data(
            manifest.key, manifest.entry.removeprefix("./")
        )
    except KeyError:
        user_entry = None
    if user_entry is not None and user_entry.is_file():
        return _status_for(
            manifest, ok=True, source=_source_for_tier(tier), path=user_entry
        )

    # 3. tools/<TOOL>/
    entry = manifest.entry_path
    if entry.is_file():
        return _status_for(manifest, ok=True, source="vendored", path=entry)

    # 4. PATH (only where path_executables lists a sensible binary)
    for exe in manifest.path_executables:
        hit = shutil.which(exe)
        if hit:
            return _status_for(manifest, ok=True, source="path", path=Path(str(hit)))

    # 5. missing
    return _status_for(manifest, ok=False, source="missing", path=None)


def contract_env(**overrides: str) -> dict[str, str]:
    """Environment for invoking a tool contract ``run`` as a child process.

    Ruling 1: a bash runner cannot discover the *caller's* ``sys.executable``,
    so every in-process invocation pins ``PHASEPRED_PYTHON`` to the interpreter
    running phasepred. ``python-library`` (LocalCIDER) and ``python`` (PScore)
    contracts then probe the same environment predict imports from, rather than
    a random ``python3`` on ``PATH``. Other contracts ignore the variable.
    Overrides are merged last so callers can still adjust a single value.
    """
    return {**os.environ, "PHASEPRED_PYTHON": sys.executable, **overrides}


def _probe_usable(entry: Path) -> bool:
    """Ask a package ``run`` script ``--check`` whether its entity is in place.

    ``PHASEPRED_PYTHON`` is pinned to the interpreter running phasepred so
    python-library contracts (LocalCIDER) probe the same env that predict
    imports from, rather than a random ``python3`` on PATH.
    """
    try:
        result = subprocess.run(
            [str(entry), "--check"],
            capture_output=True,
            text=True,
            timeout=60,
            env=contract_env(),
        )
        return result.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


# ---------------------------------------------------------------------------
# Public resolver API
# ---------------------------------------------------------------------------


def find_tool_entry(key: str) -> Path:
    """Resolve the executable to run for a tool (raises when missing)."""
    status = resolve_status(key)
    if not status.ok or status.path is None:
        raise PhaSePredToolNotFound(f"{_display(key)} not found. {status.install_hint}")
    return status.path


def resolve_status(key: str) -> ToolStatus:
    """Resolve a tool contract by key (no runtime-usability probe)."""
    manifest = _find_manifest(key)
    return _resolve(manifest)


def find_tool_entry_or_none(key: str) -> Path | None:
    """Like :func:`find_tool_entry` but returns ``None`` instead of raising."""
    try:
        return find_tool_entry(key)
    except PhaSePredToolNotFound:
        return None


# ---------------------------------------------------------------------------
# Diagnostic helpers used by `phasepred check-tools` and predict() preflight
# ---------------------------------------------------------------------------


def check_all_tools(tools_root: Path | None = None) -> list[ToolStatus]:
    """Return exactly one status row per feature-column source (9 rows).

    Six manifest contracts, the synthetic catGRANULE package row, and the
    synthetic DeepPhase and PhosphoSitePlus data rows; together they cover all
    10 feature columns (LocalCIDER supplies two). ``ok`` reflects *runtime
    usability*: the resolution must succeed and — for
    ``env``/``user``/``repo``/``vendored`` sources, where the ``run`` script is
    ours — the package's ``./run --check`` probe must pass.
    """
    statuses: list[ToolStatus] = []
    for manifest in load_manifests(tools_root):
        status = _resolve(manifest)
        if status.ok and status.path is not None and status.source != "path":
            if not _probe_usable(status.path):
                status = replace(status, ok=False)
        statuses.append(status)

    statuses.append(_catgranule_status())
    statuses.append(_deepphase_status())
    statuses.append(_phosphosite_status())
    return statuses


def _catgranule_status() -> ToolStatus:
    try:
        import catgranule as _pkg
        from catgranule import __version__ as catgranule_version
        from catgranule.scoring import score_sequence as _s  # noqa: F401

        return ToolStatus(
            name="catGRANULE",
            ok=True,
            source="package",
            path=Path(str(_pkg.__file__)).parent,
            required_for_predict=True,
            license="MIT",
            redistributable=True,
            kind="python",
            install_hint=(
                f"catgranule package v{catgranule_version} "
                "(pyproject dependency; weights artifact legacy-unaudited)."
            ),
            feature_columns=("catGRANULE",),
        )
    except Exception as e:  # pragma: no cover
        return ToolStatus(
            name="catGRANULE",
            ok=False,
            source="missing",
            path=None,
            required_for_predict=True,
            license="MIT",
            redistributable=True,
            kind="python",
            install_hint=f"catgranule package import failed: {e}",
            feature_columns=("catGRANULE",),
        )


def _deepphase_status() -> ToolStatus:
    """Synthetic row for the DeepPhase feature column.

    Preference: the wheel-bundled derived TSV (``source=package``) > the
    optional external xlsx resolved through the user-data cascade (Ruling 9:
    ``source`` is ``"user"`` only for a user-controlled root and ``"repo"`` for
    the checkout fallback). Not marked ``required_for_predict``: DeepPhase is a
    human-mode-only column, and the preflight gate is not mode-aware.
    """
    bundled = bundled_path("data", "deepphase", "deepphase_scores.tsv")
    if bundled.is_file():
        return ToolStatus(
            name="DeepPhase",
            ok=True,
            source="package",
            path=bundled,
            required_for_predict=False,
            license="MIT",
            redistributable=True,
            kind="data",
            install_hint="Bundled with the phasepred wheel (deepphase_scores.tsv).",
            feature_columns=("DeepPhase",),
        )
    xlsx, tier = resolve_optional_data("deepphase", "extracted", "tableS3.xlsx")
    if xlsx.is_file():
        return ToolStatus(
            name="DeepPhase",
            ok=True,
            source=_source_for_tier(tier),
            path=xlsx,
            required_for_predict=False,
            license="MIT",
            redistributable=True,
            kind="data",
            install_hint="Resolved from the optional external tableS3.xlsx.",
            feature_columns=("DeepPhase",),
        )
    return ToolStatus(
        name="DeepPhase",
        ok=False,
        source="missing",
        path=None,
        required_for_predict=False,
        license="MIT",
        redistributable=True,
        kind="data",
        install_hint=(
            "DeepPhase score table not found. Reinstall the phasepred wheel "
            "(bundled data/deepphase/deepphase_scores.tsv) or place tableS3.xlsx "
            f"under {user_data_root()}/deepphase/extracted/."
        ),
        feature_columns=("DeepPhase",),
    )


def _phosphosite_status() -> ToolStatus:
    """Synthetic row for the PhosphoSitePlus (Phos freq) feature column.

    Registration-gated data: never redistributed, never marked
    ``required_for_predict`` (human-mode-only column; the preflight gate is not
    mode-aware). Ruling 9: a checkout-tier hit reports ``source="repo"`` rather
    than lying with ``"user"``. The missing hint names the user data directory
    and installer.
    """
    path, tier = resolve_optional_data(
        "phosphositeplus", "Phosphorylation_site_dataset.gz"
    )
    if path.is_file():
        return ToolStatus(
            name="PhosphoSitePlus",
            ok=True,
            source=_source_for_tier(tier),
            path=path,
            required_for_predict=False,
            license="CC-BY-NC-SA-3.0",
            redistributable=False,
            kind="data",
            install_hint="Registration-gated data; do not redistribute.",
            feature_columns=("Phos freq",),
        )
    return ToolStatus(
        name="PhosphoSitePlus",
        ok=False,
        source="missing",
        path=None,
        required_for_predict=False,
        license="CC-BY-NC-SA-3.0",
        redistributable=False,
        kind="data",
        install_hint=(
            "PhosphoSitePlus dataset not found. Install it into the user data "
            f"directory: bash tools/PhosphoSitePlus/install.sh (writes to "
            f"{user_data_root()}/phosphositeplus/; override with "
            "PHASEPRED_DATA_ROOT)."
        ),
        feature_columns=("Phos freq",),
    )


def statuses_to_json(statuses: list[ToolStatus]) -> str:
    """Render status rows as a JSON array (Ruling 4: stdout stays pure JSON).

    Each object carries exactly the nine fields ``tool, status, required,
    source, path, license, redistributable, feature_columns, hint``. ``path``
    is a string or ``None`` (JSON cannot carry a ``Path``). ``source`` values
    are canonicalised through :func:`_render_source`, the same table
    :func:`format_status_table` consults. ``indent=2`` and a trailing newline.
    """
    payload = [
        {
            "tool": s.name,
            "status": "OK" if s.ok else "MISSING",
            "required": s.required_for_predict,
            "source": _render_source(s.source),
            "path": str(s.path) if s.path is not None else None,
            "license": s.license,
            "redistributable": s.redistributable,
            "feature_columns": list(s.feature_columns),
            "hint": s.install_hint,
        }
        for s in statuses
    ]
    return json.dumps(payload, indent=2) + "\n"


def assert_predict_requirements() -> list[ToolStatus]:
    """Raise PhaSePredToolNotFound if any predict-required tool is MISSING.

    Returns the full status list so callers can render it.
    """
    statuses = check_all_tools()
    missing_required = [s for s in statuses if not s.ok and s.required_for_predict]
    if missing_required:
        msg_lines = [
            "phasepred predict aborted — required external tool(s) missing:",
            "",
        ]
        for s in missing_required:
            msg_lines.append(f"  - {s.name}: {s.install_hint}")
        msg_lines.append("")
        msg_lines.append("Run `phasepred check-tools` for the full status table.")
        raise PhaSePredToolNotFound("\n".join(msg_lines))
    return statuses


def format_status_table(statuses: list[ToolStatus]) -> str:
    """Render the status list as a fixed-width table for terminal display.

    The ``SOURCE`` column renders the same :data:`SOURCE_VOCABULARY` values as
    :func:`statuses_to_json` (``env`` as ``env:<VAR>``), canonicalised through
    the shared :func:`_render_source` table; a MISSING row shows ``—`` instead,
    since its source is always ``missing``. The ``PATH`` column also shows ``—``
    for a MISSING row even when a resolution path exists, because that path is
    only the vendored ``run`` wrapper — the runtime behind it is not usable — so
    printing it beside MISSING would read as if the tool were installed.

    After the summary line, every MISSING row's ``install_hint`` is printed in a
    trailing ``Install hints`` block: the remedy already exists in the data model
    and in ``--json`` (under ``hint``), and dropping it from the table left a
    reader who only saw ``MISSING`` with no way to fix it. When any of those
    hints names a ``tools/<TOOL>/install.sh`` installer, a note states that
    ``tools/`` ships only in the GitHub repository, not in the wheel or the
    sdist, so a PyPI-only user knows the installer needs a checkout.
    """
    rows = [
        (
            "TOOL",
            "STATUS",
            "REQUIRED",
            "SOURCE",
            "PATH",
            "LICENSE",
            "REDISTRIBUTABLE",
        )
    ]
    for s in statuses:
        rows.append((
            s.name,
            "OK" if s.ok else "MISSING",
            "yes" if s.required_for_predict else "no",
            _render_source(s.source) if s.ok else "—",
            str(s.path) if s.ok and s.path else "—",
            s.license,
            "yes" if s.redistributable else "no",
        ))
    widths = [max(len(r[i]) for r in rows) for i in range(7)]
    out = []
    for idx, r in enumerate(rows):
        out.append("  ".join(r[i].ljust(widths[i]) for i in range(7)).rstrip())
        if idx == 0:
            out.append("  ".join("-" * w for w in widths))
    n_ok = sum(1 for s in statuses if s.ok)
    n_missing = sum(1 for s in statuses if not s.ok)
    n_missing_required = sum(1 for s in statuses if not s.ok and s.required_for_predict)
    out.append("")
    out.append(f"{n_ok} OK  ·  {n_missing} MISSING  ·  {n_missing_required} required-but-missing")
    missing = [s for s in statuses if not s.ok]
    if missing:
        out.append("")
        out.append("Install hints:")
        for s in missing:
            out.append(f"  {s.name}: {s.install_hint}")
        if any("tools/" in s.install_hint for s in missing):
            out.append("")
            out.append(
                "Note: `tools/` ships only in the GitHub repository (a clone), "
                "not in the wheel"
            )
            out.append(
                "or the sdist, so a hint above naming `tools/<TOOL>/install.sh` "
                "requires a checkout."
            )
    return "\n".join(out)
