set -euo pipefail

python scripts/build_FigureS4C_threshold_sweep.py
Rscript scripts/plot_FigureS4C_threshold_sweep.R

echo "threshold_sweep_metrics.tsv: $(realpath results/FigureS4C_threshold_sweep/tables/threshold_sweep_metrics.tsv)"
echo "threshold_optima.tsv: $(realpath results/FigureS4C_threshold_sweep/tables/threshold_optima.tsv)"
echo "FigureS4C_threshold_sweep.pdf: $(realpath results/FigureS4C_threshold_sweep/plots/FigureS4C_threshold_sweep.pdf)"

if [ -f results/FigureS4/plots/FigureS4.pdf ]; then
  echo "FigureS4.pdf: $(realpath results/FigureS4/plots/FigureS4.pdf)"
fi
