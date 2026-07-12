"""
Exact single-point human legacy reproduction for Phase 1.6.

This runner intentionally targets one old GAT result point only:
- organism: human
- feature subset and threshold come from config
- metric source: legacy result tables
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
    parser = argparse.ArgumentParser(description="Run exact human legacy reproduction")
    parser.add_argument("--config", required=True, type=str)
    return parser.parse_args()


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def write_yaml(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def main():
    args = parse_args()
    config = load_yaml(args.config)

    output_root = os.path.join(config["paths"]["output_root"], config["run"]["name"])
    os.makedirs(output_root, exist_ok=True)

    old_results = pd.read_csv(config["target"]["source_file"])
    target_rows = old_results[
        (old_results["organism"] == config["target"]["species"])
        & (old_results["name"] == config["target"]["row_name"])
        & (old_results["expression"].astype(bool) == False)
        & (old_results["orthologs"].astype(bool) == True)
        & (old_results["sublocs"].astype(bool) == True)
        & (old_results["string_thr"].astype(int) == int(config["target"]["string_threshold"]))
    ].copy()
    if target_rows.empty:
        raise RuntimeError("Target row not found in old results table")
    target_row = target_rows.iloc[0]

    trace_df = pd.DataFrame(
        [
            {
                "target_metric": config["target"]["target_metric"],
                "source_file": config["target"]["source_file"],
                "row_name": config["target"]["row_name"],
                "feature_combination": config["target"]["feature_combination"],
                "aggregation_method": config["target"]["aggregation_method"],
                "threshold_selection": config["target"]["threshold_selection"],
                "notes": "Exact target row chosen from legacy GAT_results-string.csv",
            }
        ]
    )
    trace_df.to_csv(os.path.join(output_root, "old_target_trace.tsv"), sep="\t", index=False)

    feature_subset = config["target"]["feature_combination"]
    old_vs_current_rows = [
        [
            {
                "step": "feature_subset",
                "old_pipeline": feature_subset,
                "current_replay": "",
                "missing_component": "exact subset selection",
                "impact": "high",
            },
            {
                "step": "string_threshold",
                "old_pipeline": str(config["target"]["string_threshold"]),
                "current_replay": "",
                "missing_component": "old threshold point",
                "impact": "high",
            },
            {
                "step": "n_runs",
                "old_pipeline": str(config["execution"]["n_runs"]),
                "current_replay": "",
                "missing_component": "aggregation level must match target",
                "impact": "medium",
            },
            {
                "step": "epochs",
                "old_pipeline": str(config["base_config"]["train"]["epochs"]),
                "current_replay": "",
                "missing_component": "old training budget",
                "impact": "high",
            },
            {
                "step": "checkpoint_selection",
                "old_pipeline": "final_model_state" if not bool(config["base_config"]["train"].get("reload_best_state", True)) else "best_val_state",
                "current_replay": "",
                "missing_component": "test-time state selection",
                "impact": "medium",
            },
            {
                "step": "summary_provenance",
                "old_pipeline": config["target"]["source_file"],
                "current_replay": "",
                "missing_component": "",
                "impact": "low",
            },
        ]
    ]
    replay_ref = config.get("current_replay_reference", {})
    if replay_ref.get("metrics_path"):
        previous_metrics = pd.read_csv(replay_ref["metrics_path"], sep="\t").iloc[0]
        old_vs_current_rows[0]["current_replay"] = replay_ref.get("feature_combination", "")
        old_vs_current_rows[1]["current_replay"] = str(replay_ref.get("string_threshold", ""))
        old_vs_current_rows[2]["current_replay"] = str(replay_ref.get("n_runs", ""))
        old_vs_current_rows[3]["current_replay"] = str(replay_ref.get("epochs", ""))
        old_vs_current_rows[4]["current_replay"] = replay_ref.get("checkpoint_selection", "")
        old_vs_current_rows[5]["current_replay"] = replay_ref.get("description", "direct replay metric")
    else:
        previous_metrics = None
    old_vs_current = pd.DataFrame(old_vs_current_rows)
    old_vs_current.to_csv(os.path.join(output_root, "old_pipeline_vs_current.tsv"), sep="\t", index=False)

    threshold_df = old_results[
        (old_results["organism"] == config["target"]["species"])
        & (old_results["name"] == config["target"]["row_name"])
        & (old_results["expression"].astype(bool) == False)
        & (old_results["orthologs"].astype(bool) == True)
        & (old_results["sublocs"].astype(bool) == True)
    ][["string_thr", "mean", "auc_pr", "accuracy"]].copy()
    threshold_df["selected_for_reproduction"] = threshold_df["string_thr"].astype(int) == int(config["target"]["string_threshold"])
    threshold_df.to_csv(os.path.join(output_root, "threshold_sweep.tsv"), sep="\t", index=False)

    run_metrics = []
    n_runs = int(config["execution"]["n_runs"])
    for run_idx in range(n_runs):
        run_seed = run_idx if config["execution"]["seed_strategy"] == "legacy_run_index" else run_idx
        run_name = os.path.join(config["run"]["name"], "runs", "run_{}".format(run_idx))
        run_config = copy.deepcopy(config["base_config"])
        run_config["run"]["name"] = run_name
        run_config["legacy"]["seed"] = run_seed
        build_dataset(run_config)
        result = train_from_config(run_config)
        row = dict(result["metrics"])
        row["run_id"] = "run_{}".format(run_idx)
        row["seed_used"] = run_seed
        run_metrics.append(row)

    all_runs = pd.DataFrame(run_metrics)
    all_runs.to_csv(os.path.join(output_root, "all_runs.tsv"), sep="\t", index=False)

    subprocess.check_call(
        [
            "python",
            "src/eval/aggregate_epgat_legacy_runs.py",
            "--run-root",
            output_root,
        ]
    )

    aggregated = pd.read_csv(os.path.join(output_root, "aggregated_metrics.tsv"), sep="\t")
    aggregated.to_csv(os.path.join(output_root, "aggregated.tsv"), sep="\t", index=False)
    agg_row = aggregated.iloc[0]

    best_selected = pd.DataFrame(
        [
            {
                "metric": "auroc_mean",
                "selected_threshold": int(config["target"]["string_threshold"]),
                "row_name": config["target"]["row_name"],
                "old_value": float(target_row["mean"]),
                "reproduced_value": float(agg_row["auroc_mean"]),
            }
        ]
    )
    best_selected.to_csv(os.path.join(output_root, "best_selected_result.tsv"), sep="\t", index=False)

    exact_vs_old = pd.DataFrame(
        [
            {
                "metric": "auroc",
                "old_value": float(target_row["mean"]),
                "reproduced_value": float(agg_row["auroc_mean"]),
                "delta": float(agg_row["auroc_mean"]) - float(target_row["mean"]),
                "match_status": "match" if abs(float(agg_row["auroc_mean"]) - float(target_row["mean"])) <= 0.01 else ("close" if abs(float(agg_row["auroc_mean"]) - float(target_row["mean"])) <= 0.03 else "mismatch"),
            },
            {
                "metric": "auprc",
                "old_value": float(target_row["auc_pr"]),
                "reproduced_value": float(agg_row["auprc_mean"]),
                "delta": float(agg_row["auprc_mean"]) - float(target_row["auc_pr"]),
                "match_status": "match" if abs(float(agg_row["auprc_mean"]) - float(target_row["auc_pr"])) <= 0.02 else ("close" if abs(float(agg_row["auprc_mean"]) - float(target_row["auc_pr"])) <= 0.05 else "mismatch"),
            },
            {
                "metric": "accuracy",
                "old_value": float(target_row["accuracy"]),
                "reproduced_value": float(agg_row["accuracy_mean"]),
                "delta": float(agg_row["accuracy_mean"]) - float(target_row["accuracy"]),
                "match_status": "match" if abs(float(agg_row["accuracy_mean"]) - float(target_row["accuracy"])) <= 0.02 else ("close" if abs(float(agg_row["accuracy_mean"]) - float(target_row["accuracy"])) <= 0.05 else "mismatch"),
            },
        ]
    )
    exact_vs_old.to_csv(os.path.join(output_root, "exact_vs_old.tsv"), sep="\t", index=False)

    write_yaml(os.path.join(output_root, "run_config_frozen.yaml"), config)

    lines = [
        "# Exact Reproduction Summary",
        "",
        "- target row: {}".format(config["target"]["row_name"]),
        "- target threshold: {}".format(config["target"]["string_threshold"]),
        "- old AUROC: {:.4f}".format(float(target_row["mean"])),
        "- reproduced AUROC: {:.4f}".format(float(agg_row["auroc_mean"])),
        "- old AUPRC: {:.4f}".format(float(target_row["auc_pr"])),
        "- reproduced AUPRC: {:.4f}".format(float(agg_row["auprc_mean"])),
        "- old Accuracy: {:.4f}".format(float(target_row["accuracy"])),
        "- reproduced Accuracy: {:.4f}".format(float(agg_row["accuracy_mean"])),
    ]
    if previous_metrics is not None:
        lines.append("- previous replay AUROC: {:.4f}".format(float(previous_metrics["auroc"])))
    with open(os.path.join(output_root, "exact_reproduction_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


if __name__ == "__main__":
    main()
