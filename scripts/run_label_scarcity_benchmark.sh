#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/jiehuang/software/fungi/ProGATE_v2"
CORES="4"
declare -a EXTRA_ARGS=()
declare -a TARGETS=()

usage() {
  cat <<'EOF'
Usage:
  scripts/run_label_scarcity_benchmark.sh [options] [targets...]

Options:
  -c, --cores N   Set Snakemake cores. Default is 4.
  -h, --help      Show this help message.

Examples:
  scripts/run_label_scarcity_benchmark.sh --cores 16
  scripts/run_label_scarcity_benchmark.sh --cores 8 summarize_label_scarcity
EOF
}

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
mkdir -p "${PROJECT_ROOT}/.mplconfig" "${PROJECT_ROOT}/.cache" "${PROJECT_ROOT}/.tmp" "outputs/Figure2_label_scarcity" "results/Figure2_label_scarcity"

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
  --snakefile workflow/label_scarcity_benchmark.smk
  --cores "${CORES}"
  --rerun-incomplete
  --keep-going
)

if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  CMD+=("${EXTRA_ARGS[@]}")
fi

if [[ ${#TARGETS[@]} -gt 0 ]]; then
  CMD+=(-- "${TARGETS[@]}")
fi

echo "[label_scarcity] protocol=fgraminearum_newlabel (default; override via Snakemake --config if needed)"
echo "[label_scarcity] feature_setting=workflow-configured"
echo "[label_scarcity] output_root=outputs/Figure2_label_scarcity"
echo "[label_scarcity] results_root=results/Figure2_label_scarcity"

exec env XDG_CACHE_HOME="${XDG_CACHE_HOME}" TMPDIR="${TMPDIR}" "${CMD[@]}"
