#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"
export PYTHONPATH="${ROOT_DIR}:${PYTHONPATH:-}"
export MPLCONFIGDIR="${ROOT_DIR}/results/Figure3c_threshold_tuned/.matplotlib"

python src/eval/summarize_figure3c_threshold_tuned_metrics.py \
  --figure3c-root outputs/Figure3c/fgraminearum_newlabel/GraphSAGE \
  --summary-dir results/Figure3c_threshold_tuned

python src/plot/plot_figure3c_baseline_vs_gated_roc_pr.py \
  --figure3c-root outputs/Figure3c/fgraminearum_newlabel/GraphSAGE \
  --output-dir results/Figure3c_threshold_tuned/plots \
  --mode pooled

python src/plot/plot_figure3c_baseline_vs_gated_roc_pr.py \
  --figure3c-root outputs/Figure3c/fgraminearum_newlabel/GraphSAGE \
  --output-dir results/Figure3c_threshold_tuned/plots \
  --mode mean_seed
