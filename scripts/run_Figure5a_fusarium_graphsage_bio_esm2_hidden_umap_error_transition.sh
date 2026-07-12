#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${ROOT_DIR}/.mplconfig}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/codex_cache}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp/numba_cache}"
export TMPDIR="${TMPDIR:-/tmp}"

echo "[Figure5a] hidden_umap_error_transition"
echo "[Figure5a] workflow=workflow/Figure5a_hidden_umap_error_transition.smk"
echo "[Figure5a] results_root=results/Figure5"

has_cores_arg=false
for arg in "$@"; do
  case "$arg" in
    -j|--cores|--cores=*|-j*)
      has_cores_arg=true
      break
      ;;
  esac
done

if [[ "${has_cores_arg}" == true ]]; then
  snakemake --snakefile "workflow/Figure5a_hidden_umap_error_transition.smk" figure5a "$@"
else
  snakemake --snakefile "workflow/Figure5a_hidden_umap_error_transition.smk" --cores 1 figure5a "$@"
fi
