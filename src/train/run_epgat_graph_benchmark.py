"""
Run the archived 4-species x 4-model replay benchmark.
"""

import argparse
import os
import subprocess


CONFIGS = [
    "configs/epgat_graph_benchmark_human.yaml",
    "configs/epgat_graph_benchmark_celegans.yaml",
    "configs/epgat_graph_benchmark_scerevisiae.yaml",
    "configs/epgat_graph_benchmark_fgraminearum.yaml",
]
MODELS = ["gat", "gcn", "gin", "sage"]


def main():
    parser = argparse.ArgumentParser(description="Run graph benchmark")
    parser.add_argument("--output-root", default="outputs/epgat_graph_benchmark", type=str)
    parser.parse_args()

    for config in CONFIGS:
        species = os.path.basename(config).replace("configs/epgat_graph_benchmark_", "").replace(".yaml", "")
        for model in MODELS:
            metrics_path = os.path.join("outputs/epgat_graph_benchmark", species, model, "metrics.tsv")
            if os.path.exists(metrics_path):
                continue
            subprocess.check_call(["python", "src/train/train_epgat_graph_models.py", "--config", config, "--model", model])


if __name__ == "__main__":
    main()
