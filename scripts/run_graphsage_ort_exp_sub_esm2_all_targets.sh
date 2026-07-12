#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${REPO_ROOT}"

CONFIG="configs/frozen_protocol.yaml"
SNAKEFILE="workflow/frozen_protocol_benchmark.smk"
OUTPUT_ROOT="outputs/frozen_protocol_benchmark_v2"
SUMMARY_DIR="results/frozen_protocol/graphsage_ort_exp_sub_esm2_all_targets"
SUMMARY_PYTHON_MODULE="src.eval.summarize_graphsage_ort_exp_sub_esm2_all_targets"
MODEL="GraphSAGE_ORT_EXP_SUB_ESM2"
FEATURE="ORT_EXP_SUB_ESM2"
PROTOCOLS=(
  "human"
  "celegans"
  "scerevisiae"
  "dmelanogaster"
  "fgraminearum_oldlabel"
  "fgraminearum_newlabel"
)

DEFAULT_CORES="4"
CORES="${DEFAULT_CORES}"
DRY_RUN="false"
SEED_FILTER=""

usage() {
  cat <<'EOF'
Usage:
  scripts/run_graphsage_ort_exp_sub_esm2_all_targets.sh [options] [cores]

Options:
  --cores N      Set Snakemake cores. If omitted, default is 4.
  --seed N       Run only one seed as a smoke test. Default is all five frozen seeds.
  --dry-run      Run Snakemake in dry-run mode and add --nolock.
  -h, --help     Show this help message.

Examples:
  scripts/run_graphsage_ort_exp_sub_esm2_all_targets.sh --dry-run
  scripts/run_graphsage_ort_exp_sub_esm2_all_targets.sh --seed 1029 --cores 1
  scripts/run_graphsage_ort_exp_sub_esm2_all_targets.sh 8
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cores)
      [[ $# -ge 2 ]] || { echo "Missing value for --cores" >&2; exit 1; }
      CORES="$2"
      shift 2
      ;;
    --seed)
      [[ $# -ge 2 ]] || { echo "Missing value for --seed" >&2; exit 1; }
      SEED_FILTER="$2"
      shift 2
      ;;
    --dry-run)
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
      if [[ "${CORES}" != "${DEFAULT_CORES}" ]]; then
        echo "Positional cores cannot be combined with --cores" >&2
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

if [[ -n "${SEED_FILTER}" ]] && ! [[ "${SEED_FILTER}" =~ ^[0-9]+$ ]]; then
  echo "Seed must be an integer, got: ${SEED_FILTER}" >&2
  exit 1
fi

if [[ -z "${HOME:-}" || ! -w "${HOME}" ]]; then
  export HOME="${REPO_ROOT}"
fi

export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${REPO_ROOT}/.cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${REPO_ROOT}/.mplconfig}"
export TMPDIR="${TMPDIR:-/tmp}"
mkdir -p "${XDG_CACHE_HOME}" "${MPLCONFIGDIR}"

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base 2>/dev/null || true)"
  if [[ -n "${CONDA_BASE}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    if conda env list | awk '{print $1}' | grep -qx "progate"; then
      conda activate progate || true
    else
      echo "[graphsage_ort_exp_sub_esm2_all_targets] conda env 'progate' not found; continuing without activation" >&2
    fi
  fi
fi

export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export MPLBACKEND="Agg"

if [[ -n "${SEED_FILTER}" ]]; then
  SEEDS=("${SEED_FILTER}")
else
  SEEDS=("1029" "1030" "1031" "1032" "1033")
fi

TARGETS=()
for protocol in "${PROTOCOLS[@]}"; do
  for seed in "${SEEDS[@]}"; do
    TARGETS+=("${OUTPUT_ROOT}/${protocol}/${MODEL}/${FEATURE}/run_${seed}/metrics.tsv")
  done
done

echo "[graphsage_ort_exp_sub_esm2_all_targets] repo_root=${REPO_ROOT}"
echo "[graphsage_ort_exp_sub_esm2_all_targets] model=${MODEL}"
echo "[graphsage_ort_exp_sub_esm2_all_targets] feature=${FEATURE}"
echo "[graphsage_ort_exp_sub_esm2_all_targets] protocols=${PROTOCOLS[*]}"
echo "[graphsage_ort_exp_sub_esm2_all_targets] seeds=${SEEDS[*]}"
echo "[graphsage_ort_exp_sub_esm2_all_targets] targets=${#TARGETS[@]}"
echo "[graphsage_ort_exp_sub_esm2_all_targets] summary_dir=${SUMMARY_DIR}"
printf '  - %s\n' "${TARGETS[@]}"

SNAKEMAKE_ARGS=(
  -s "${SNAKEFILE}"
  --configfile "${CONFIG}"
  --cores "${CORES}"
  --rerun-incomplete
  --keep-going
  --rerun-triggers=mtime
)

if [[ "${DRY_RUN}" == "true" ]]; then
  SNAKEMAKE_ARGS+=(--dry-run --nolock)
fi

if command -v snakemake >/dev/null 2>&1; then
  SNAKEMAKE_CMD=(snakemake)
else
  SNAKEMAKE_CMD=(python -m snakemake)
fi

run_snakemake() {
  "${SNAKEMAKE_CMD[@]}" "${SNAKEMAKE_ARGS[@]}" "${TARGETS[@]}"
}

unlock_snakemake() {
  XDG_CACHE_HOME="${XDG_CACHE_HOME}" "${SNAKEMAKE_CMD[@]}" \
    -s "${SNAKEFILE}" \
    --configfile "${CONFIG}" \
    --unlock
}

if ! run_snakemake; then
  if [[ "${DRY_RUN}" == "true" ]]; then
    exit 1
  fi
  echo "[graphsage_ort_exp_sub_esm2_all_targets] snakemake failed; attempting one automatic unlock + retry" >&2
  unlock_snakemake
  run_snakemake
fi

if [[ "${DRY_RUN}" == "true" ]]; then
  echo "[graphsage_ort_exp_sub_esm2_all_targets] dry-run mode: skipped summary generation"
  echo "[graphsage_ort_exp_sub_esm2_all_targets] summary outputs would be written under ${SUMMARY_DIR}"
  exit 0
fi

mkdir -p "${SUMMARY_DIR}"
python -m "${SUMMARY_PYTHON_MODULE}" \
  --output-root "${OUTPUT_ROOT}" \
  --summary-dir "${SUMMARY_DIR}"

echo "[graphsage_ort_exp_sub_esm2_all_targets] summary written to ${SUMMARY_DIR}"
