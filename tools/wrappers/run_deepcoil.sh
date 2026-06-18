#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
env_prefix="${DEEPCOIL_ENV_PREFIX:-$repo_root/.external_envs/deepcoil}"

if [[ ! -x "$env_prefix/bin/deepcoil" ]]; then
  cat >&2 <<EOF
DeepCoil executable not found: $env_prefix/bin/deepcoil

Create it with:
  tools/install/install_deepcoil_env.sh
EOF
  exit 1
fi

exec "$env_prefix/bin/deepcoil" "$@"
