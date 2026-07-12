"""
Train the slim Phase 2 extended dataset using the existing legacy GAT.
"""

import argparse
import copy

import yaml

from src.train.train_epgat_legacy import train_from_config


def parse_args():
    parser = argparse.ArgumentParser(description="Train slim extended EPGAT run")
    parser.add_argument("--config", required=True, type=str)
    parser.add_argument("--feature-mode", required=True, choices=["baseline", "baseline_plus_esm2", "baseline_plus_prott5"])
    return parser.parse_args()


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def main():
    args = parse_args()
    config = load_yaml(args.config)
    train_config = copy.deepcopy(config["train_config"])
    train_config["run"]["name"] = "{}__{}".format(config["run"]["name"], args.feature_mode)
    train_config["paths"]["output_root"] = config["paths"]["output_root"]
    result = train_from_config(train_config)
    print("Slim extended training complete:", result["run_out"])
    print(result["metrics"])


if __name__ == "__main__":
    main()
