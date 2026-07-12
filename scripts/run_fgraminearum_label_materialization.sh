#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/jiehuang/software/fungi/ProGATE_v2"
cd "${PROJECT_ROOT}"

CONFIG="configs/fgraminearum_label_materialization.yaml"
CORES="${1:-4}"
if [[ $# -gt 0 ]]; then
  shift
fi

mkdir -p "${PROJECT_ROOT}/.mplconfig" "${PROJECT_ROOT}/.cache"

if [ -f "/home/jiehuang/anaconda3/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  source "/home/jiehuang/anaconda3/etc/profile.d/conda.sh"
  conda activate EPGAT
fi

export MPLBACKEND="Agg"
export MPLCONFIGDIR="${PROJECT_ROOT}/.mplconfig"
export XDG_CACHE_HOME="${PROJECT_ROOT}/.cache"

exec /home/jiehuang/anaconda3/bin/snakemake \
  -s workflow/fgraminearum_label_materialization.smk \
  --cores "${CORES}" \
  --rerun-incomplete \
  --keep-going \
  --configfile "${CONFIG}" \
  "$@"
