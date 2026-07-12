#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/jiehuang/software/fungi/ProGATE_v2"
cd "${PROJECT_ROOT}"

CONFIG="configs/frozen_protocol.yaml"
DEFAULT_CORES="4"
CORES="${DEFAULT_CORES}"

usage() {
  cat <<'EOF'
Usage:
  scripts/run_Figure1_frozen_protocol_benchmark.sh [options] [targets...]

Options:
  -c, --cores N  Set Snakemake cores. If omitted, default is 4.
  -h, --help     Show this help message.

Notes:
  Remaining options starting with '-' are passed through to Snakemake.
  Positional non-option arguments are treated as explicit Snakemake targets.

Examples:
  scripts/run_Figure1_frozen_protocol_benchmark.sh -c 64
  scripts/run_Figure1_frozen_protocol_benchmark.sh --cores 8 --dry-run
  scripts/run_Figure1_frozen_protocol_benchmark.sh results/Figure1/summary/final_summary.tsv
EOF
}

declare -a SNAKEMAKE_ARGS=()
declare -a TARGETS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--cores)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 1; }
      CORES="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while [[ $# -gt 0 ]]; do
        TARGETS+=("$1")
        shift
      done
      ;;
    -*)
      SNAKEMAKE_ARGS+=("$1")
      shift
      if [[ $# -gt 0 && "$1" != -* ]]; then
        case "${SNAKEMAKE_ARGS[-1]}" in
          --batch|--config|--configfile|--default-resources|--forceuse-threads|--groups|--group-components|--jobs|--list-changes|--list-input-changes|--list-params-changes|--list-untracked|--max-inventory-time|--max-jobs-per-second|--max-status-checks-per-second|--resources|--scheduler|--set-threads|--set-resources|--set-scatter|--snakefile|--target-jobs|--until|--omit-from|--rerun-triggers|--allowed-rules|--local-cores|--mode|--wms-monitor-arg|--target-files-omit-workdir-adjustment)
            SNAKEMAKE_ARGS+=("$1")
            shift
            ;;
        esac
      fi
      ;;
    *)
      if [[ "${CORES}" == "${DEFAULT_CORES}" && "$1" =~ ^[0-9]+$ && ${#TARGETS[@]} -eq 0 ]]; then
        CORES="$1"
      else
        TARGETS+=("$1")
      fi
      shift
      ;;
  esac
done

if ! [[ "${CORES}" =~ ^[0-9]+$ ]] || [[ "${CORES}" -lt 1 ]]; then
  echo "Cores must be a positive integer, got: ${CORES}" >&2
  exit 1
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

echo "[Figure1 benchmark] running workflow/frozen_protocol_benchmark.smk with config ${CONFIG}" >&2

BOOTSTRAP_TARGETS=(
  results/Figure1/labels/label_protocol_summary.md
  results/Figure1/labels/fgraminearum_label_protocol_summary.md
  results/Figure1/splits/split_protocol_summary.md
  results/Figure1/splits/fgraminearum_split_protocol_summary.md
)

SKIP_BOOTSTRAP=false
if [[ ${#SNAKEMAKE_ARGS[@]} -gt 0 ]]; then
  for arg in "${SNAKEMAKE_ARGS[@]}"; do
    case "${arg}" in
      --unlock|--cleanup-metadata|--list-*|--summary|--detailed-summary|--lint|--dag|--rulegraph|--filegraph)
        SKIP_BOOTSTRAP=true
        ;;
    esac
  done
fi

if [[ "${SKIP_BOOTSTRAP}" == false ]]; then
  /home/jiehuang/anaconda3/bin/snakemake \
    -s workflow/frozen_protocol_benchmark.smk \
    --cores 1 \
    --configfile "${CONFIG}" \
    --rerun-incomplete \
    --keep-going \
    -- "${BOOTSTRAP_TARGETS[@]}"
fi

CMD=(
  /home/jiehuang/anaconda3/bin/snakemake
  -s workflow/frozen_protocol_benchmark.smk
  --cores "${CORES}"
  --rerun-incomplete
  --keep-going
  --configfile "${CONFIG}"
)

if [[ ${#SNAKEMAKE_ARGS[@]} -gt 0 ]]; then
  CMD+=("${SNAKEMAKE_ARGS[@]}")
fi

if [[ ${#TARGETS[@]} -gt 0 ]]; then
  CMD+=(-- "${TARGETS[@]}")
fi

exec "${CMD[@]}"
