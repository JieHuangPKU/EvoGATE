#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${ROOT_DIR}/.mplconfig}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/codex_cache}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp/numba_cache}"
export TMPDIR="${TMPDIR:-/tmp}"

echo "[Figure5c] input_vs_hidden_compare"
echo "[Figure5c] workflow=workflow/Figure5c_input_vs_hidden_compare.smk"
echo "[Figure5c] results_root=results/Figure5"

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
  snakemake --snakefile "workflow/Figure5c_input_vs_hidden_compare.smk" figure5c "$@"
else
  snakemake --snakefile "workflow/Figure5c_input_vs_hidden_compare.smk" --cores 1 figure5c "$@"
fi
