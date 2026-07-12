#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

export MPLBACKEND="${MPLBACKEND:-Agg}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${ROOT_DIR}/.mplconfig}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/codex_cache}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp/numba_cache}"
export TMPDIR="${TMPDIR:-/tmp}"

echo "[Figure5] representation_export"
echo "[Figure5] status=obsolete"
echo "[Figure5] legacy export-style representation plots are no longer part of the active Figure5 manuscript workflow"
echo "[Figure5] use ./scripts/run_Figure5_fusarium_graphsage_bio_esm2_representation_mechanism.sh -j 48 instead"
exit 1
