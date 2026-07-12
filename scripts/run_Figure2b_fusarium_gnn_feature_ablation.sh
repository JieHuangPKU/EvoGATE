#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/jiehuang/software/fungi/ProGATE_v2"
CORES="4"
DRY_RUN="false"

usage() {
  cat <<'EOF'
Usage:
  scripts/run_Figure2b_fusarium_gnn_feature_ablation.sh [cores] [-n|--dry-run]
  scripts/run_Figure2b_fusarium_gnn_feature_ablation.sh [-n|--dry-run] [cores]

Options:
  -n, --dry-run   Run Snakemake in dry-run mode.
  -h, --help      Show this help message.

Examples:
  scripts/run_Figure2b_fusarium_gnn_feature_ablation.sh 48
  scripts/run_Figure2b_fusarium_gnn_feature_ablation.sh 48 -n
  scripts/run_Figure2b_fusarium_gnn_feature_ablation.sh -n 48
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--dry-run)
      DRY_RUN="true"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "Unsupported option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      if [[ "${CORES}" != "4" ]]; then
        echo "Cores already set to ${CORES}; unexpected extra positional argument: $1" >&2
        usage >&2
        exit 1
      fi
      CORES="$1"
      shift
      ;;
  esac
done

if ! [[ "${CORES}" =~ ^[0-9]+$ ]] || [[ "${CORES}" -lt 1 ]]; then
  echo "Cores must be a positive integer, got: ${CORES}" >&2
  exit 1
fi

cd "${PROJECT_ROOT}"

mkdir -p "${PROJECT_ROOT}/.mplconfig" "${PROJECT_ROOT}/.cache" "${PROJECT_ROOT}/.tmp"

if [ -f "/home/jiehuang/anaconda3/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  source "/home/jiehuang/anaconda3/etc/profile.d/conda.sh"
  conda activate EPGAT
fi

export PYTHONPATH="${PYTHONPATH:-.}:."
export MPLBACKEND="Agg"
export MPLCONFIGDIR="${PROJECT_ROOT}/.mplconfig"
export XDG_CACHE_HOME="${PROJECT_ROOT}/.cache"
export TMPDIR="${PROJECT_ROOT}/.tmp"

SNAKEMAKE_ARGS=(
  --snakefile workflow/Figure2b_fusarium_gnn_feature_ablation.smk
  --cores "${CORES}"
  --rerun-incomplete
  --keep-going
)

if [[ "${DRY_RUN}" == "true" ]]; then
  SNAKEMAKE_ARGS+=(--dry-run --nolock)
fi

exec env XDG_CACHE_HOME="${XDG_CACHE_HOME}" TMPDIR="${TMPDIR}" /home/jiehuang/anaconda3/bin/snakemake "${SNAKEMAKE_ARGS[@]}"
