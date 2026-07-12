#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/jiehuang/software/fungi/ProGATE_v2"
CORES="${1:-4}"

cd "${PROJECT_ROOT}"

mkdir -p "${PROJECT_ROOT}/.mplconfig" "${PROJECT_ROOT}/.cache"

if [ -f "/home/jiehuang/anaconda3/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  source "/home/jiehuang/anaconda3/etc/profile.d/conda.sh"
  conda activate EPGAT
fi

export PYTHONPATH="${PYTHONPATH:-.}:."
export MPLBACKEND="Agg"
export MPLCONFIGDIR="${PROJECT_ROOT}/.mplconfig"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${PROJECT_ROOT}/.cache}"
export TMPDIR="${TMPDIR:-${PROJECT_ROOT}/.tmp}"
mkdir -p "${XDG_CACHE_HOME}" "${TMPDIR}"

exec env XDG_CACHE_HOME="${XDG_CACHE_HOME}" TMPDIR="${TMPDIR}" \
  /home/jiehuang/anaconda3/bin/snakemake \
  --snakefile workflow/Figure3cA_fusarium_graphsage_ort_exp_sub_esm2_gated_residual.smk \
  --cores "${CORES}" \
  --rerun-incomplete \
  --keep-going
