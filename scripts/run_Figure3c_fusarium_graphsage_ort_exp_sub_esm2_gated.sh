#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/jiehuang/software/fungi/ProGATE_v2"
CORES="${1:-4}"
DRY_RUN="false"
EXTRA_ARGS=()

shift_count=0
if [ $# -ge 1 ]; then
  shift_count=1
fi
if [ $# -ge 2 ]; then
  case "${2}" in
    -n|--dry-run|dry-run)
      DRY_RUN="true"
      shift_count=2
      ;;
  esac
fi

if [ $# -gt "${shift_count}" ]; then
  EXTRA_ARGS=("${@:$((${shift_count}+1))}")
fi

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

SNAKE_ARGS=(
  --snakefile workflow/Figure3c_fusarium_graphsage_ort_exp_sub_esm2_gated.smk
  --cores "${CORES}"
  --rerun-incomplete
  --keep-going
)

if [ "${DRY_RUN}" = "true" ]; then
  SNAKE_ARGS+=(--dry-run --nolock)
fi

CMD=(/home/jiehuang/anaconda3/bin/snakemake "${SNAKE_ARGS[@]}")
if [ ${#EXTRA_ARGS[@]} -gt 0 ]; then
  CMD+=("${EXTRA_ARGS[@]}")
fi

exec env XDG_CACHE_HOME="${XDG_CACHE_HOME}" TMPDIR="${TMPDIR}" "${CMD[@]}"
