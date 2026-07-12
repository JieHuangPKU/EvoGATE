"""
Run a minimal Phase 1.6 fidelity-restored legacy replay.

This wrapper rebuilds per-run datasets with old-style split semantics and
repeats training for a small set of species, then aggregates the outputs.
"""

import argparse
import copy
import os
import subprocess

import pandas as pd
import yaml

from src.data.build_epgat_legacy_dataset import build_dataset
from src.train.train_epgat_legacy import train_from_config


def parse_args():
    parser = argparse.ArgumentParser(description="Run legacy replay fidelity restoration")
    parser.add_argument("--config", required=True, type=str)
    return parser.parse_args()


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def write_yaml(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def build_config_diff(config, old_condition_rows):
    rows = []
    prev = config["previous_replay"]
    restored = config["fidelity"]
    for row in old_condition_rows:
        rows.append(
            {
                "parameter": row["parameter"],
                "old_condition": row["old_condition"],
                "previous_replay_condition": prev.get(row["parameter"], ""),
                "fidelity_restored_condition": restored.get(row["parameter"], ""),
                "status": row["status"],
                "notes": row["notes"],
            }
        )
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    config = load_yaml(args.config)
    root_dir = os.path.join(config["paths"]["output_root"], config["run"]["name"])
    os.makedirs(root_dir, exist_ok=True)

    old_condition_rows = config["old_condition_rows"]
    diff_df = build_config_diff(config, old_condition_rows)
    diff_df.to_csv(os.path.join(root_dir, "config_diff_vs_old.tsv"), sep="\t", index=False)

    run_metrics = []
    previous_metrics = pd.read_csv(config["previous"]["metrics_path"], sep="\t").iloc[0].to_dict()

    n_runs = int(config["fidelity"]["n_runs"])
    epochs = int(config["fidelity"]["epochs"])
    old_seed_base = int(config["fidelity"]["seed_base"])
    old_seed_strategy = config["fidelity"]["seed_strategy"]

    for run_idx in range(n_runs):
        if old_seed_strategy == "legacy_run_index":
            run_seed = run_idx
        else:
            run_seed = old_seed_base + run_idx

        run_name = os.path.join(config["run"]["name"], "runs", "run_{}".format(run_idx))
        run_config = copy.deepcopy(config["base_config"])
        run_config["run"]["name"] = run_name
        run_config["legacy"]["seed"] = run_seed
        run_config["legacy"]["test_fraction"] = float(config["fidelity"]["test_fraction"])
        run_config["legacy"]["val_fraction"] = float(config["fidelity"]["val_fraction_total"])
        run_config["train"]["epochs"] = epochs
        run_config["train"]["dropout"] = float(config["fidelity"]["dropout"])
        run_config["train"]["reload_best_state"] = bool(config["fidelity"]["reload_best_state"])

        build_dataset(run_config)
        result = train_from_config(run_config)
        metric_row = dict(result["metrics"])
        metric_row["run_id"] = "run_{}".format(run_idx)
        metric_row["seed_used"] = run_seed
        run_metrics.append(metric_row)

    all_metrics = pd.DataFrame(run_metrics)
    all_metrics.to_csv(os.path.join(root_dir, "all_run_metrics.tsv"), sep="\t", index=False)

    subprocess.check_call(
        [
            "python",
            "src/eval/aggregate_epgat_legacy_runs.py",
            "--run-root",
            root_dir,
        ]
    )

    aggregated = pd.read_csv(os.path.join(root_dir, "aggregated_metrics.tsv"), sep="\t").iloc[0].to_dict()

    frozen = copy.deepcopy(config)
    frozen["executed_runs"] = run_metrics
    write_yaml(os.path.join(root_dir, "run_config_frozen.yaml"), frozen)

    improvement_rows = []
    for metric in ["auroc", "auprc", "mcc", "accuracy", "f1"]:
        restored_value = aggregated.get(metric + "_mean", "")
        previous_value = previous_metrics.get(metric, "")
        if restored_value == "":
            continue
        improvement_rows.append(
            {
                "metric": metric,
                "previous_replay": previous_value,
                "fidelity_restored_replay": restored_value,
                "delta": float(restored_value) - float(previous_value),
                "notes": "mean over {} runs".format(n_runs),
            }
        )
    pd.DataFrame(improvement_rows).to_csv(
        os.path.join(root_dir, "fidelity_improvement_vs_previous.tsv"), sep="\t", index=False
    )

    compare_rows = []
    for row in config["old_reference_rows"]:
        metric = row["metric_or_condition"]
        if metric in ["auroc", "auprc", "accuracy", "f1", "mcc"]:
            restored_value = aggregated.get(metric + "_mean", "")
            gap = "" if restored_value == "" else float(restored_value) - float(row["old_summary_value"])
        else:
            restored_value = row.get("restored_replay_value_override", "")
            gap = ""
        compare_rows.append(
            {
                "metric_or_condition": metric,
                "old_summary_value": row["old_summary_value"],
                "restored_replay_value": restored_value,
                "gap": gap,
                "interpretation": row["interpretation"],
            }
        )
    pd.DataFrame(compare_rows).to_csv(
        os.path.join(root_dir, "restored_vs_old_summary.tsv"), sep="\t", index=False
    )


if __name__ == "__main__":
    main()
