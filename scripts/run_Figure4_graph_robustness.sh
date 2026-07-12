#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/jiehuang/software/fungi/ProGATE_v2"
CORES="4"
THRESHOLD_MODE="main"
declare -a EXTRA_ARGS=()
declare -a TARGETS=()

usage() {
  cat <<'EOF'
Usage:
  scripts/run_Figure4_graph_robustness.sh [options] [targets...]

Options:
  -c, --cores N   Set Snakemake cores. Default is 4.
  --quick         Optional debug mode. Runs only STRING 300 and 700, but still builds main/supplementary outputs.
  --main          Explicitly run the formal mainline thresholds 100,200,...,900 (default).
  -h, --help      Show this help message.

Notes:
  - Default mode is the formal Figure4 mainline.
  - Extra Snakemake flags such as --unlock are passed through unchanged.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--cores)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 1; }
      CORES="$2"
      shift 2
      ;;
    --quick)
      THRESHOLD_MODE="quick"
      shift
      ;;
    --main)
      THRESHOLD_MODE="main"
      shift
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
    -* )
      EXTRA_ARGS+=("$1")
      shift
      ;;
    *)
      TARGETS+=("$1")
      shift
      ;;
  esac
done

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/.mplconfig" "${PROJECT_ROOT}/.cache" "${PROJECT_ROOT}/.tmp" "outputs/Figure4" "results/Figure4"

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

CMD=(
  /home/jiehuang/anaconda3/bin/snakemake
  --snakefile workflow/Figure4_graph_robustness.smk
  --cores "${CORES}"
  --rerun-incomplete
  --keep-going
  --config "figure4_threshold_mode=${THRESHOLD_MODE}"
)

if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  CMD+=("${EXTRA_ARGS[@]}")
fi

if [[ ${#TARGETS[@]} -gt 0 ]]; then
  CMD+=(-- "${TARGETS[@]}")
fi

echo "[Figure4_graph_robustness] protocol=fgraminearum_newlabel"
echo "[Figure4_graph_robustness] model=GraphSAGE"
echo "[Figure4_graph_robustness] feature_setting=ORT_EXP_SUB_ESM2"
echo "[Figure4_graph_robustness] threshold_mode=${THRESHOLD_MODE}"
echo "[Figure4_graph_robustness] results_root=results/Figure4"
echo "[Figure4_graph_robustness] output_root=outputs/Figure4"

exec env XDG_CACHE_HOME="${XDG_CACHE_HOME}" TMPDIR="${TMPDIR}" "${CMD[@]}"
