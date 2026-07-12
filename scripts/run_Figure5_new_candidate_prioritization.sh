#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${ROOT_DIR}/.mplconfig}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/codex_cache}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp/numba_cache}"
export TMPDIR="${TMPDIR:-/tmp}"
export PYTHONPATH="${PYTHONPATH:-.}:."

OUTPUT_ROOT="results/Figure5_new_candidate_prioritization"
mkdir -p "${OUTPUT_ROOT}" "${OUTPUT_ROOT}/archive_old_hidden_umap"

python -m src.eval.build_figure5_candidate_prioritization --output-root "${OUTPUT_ROOT}" "$@"
