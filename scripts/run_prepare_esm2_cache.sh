#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${REPO_ROOT}"

if [[ -z "${HOME:-}" || ! -w "${HOME}" ]]; then
  export HOME="${REPO_ROOT}"
fi

export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${REPO_ROOT}/.cache}"
mkdir -p "${XDG_CACHE_HOME}"
export HF_HOME="${HF_HOME:-${XDG_CACHE_HOME}/huggingface}"
mkdir -p "${HF_HOME}"

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base 2>/dev/null || true)"
  if [[ -n "${CONDA_BASE}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    if conda env list | awk '{print $1}' | grep -qx "progate"; then
      conda activate progate || true
    else
      echo "[prepare_esm2_cache] conda env 'progate' not found; continuing without activation" >&2
    fi
  fi
fi

export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

if command -v snakemake >/dev/null 2>&1; then
  snakemake -s workflow/prepare_esm2_cache.smk --configfile configs/prepare_esm2_cache.yaml "$@"
else
  python -m snakemake -s workflow/prepare_esm2_cache.smk --configfile configs/prepare_esm2_cache.yaml "$@"
fi
