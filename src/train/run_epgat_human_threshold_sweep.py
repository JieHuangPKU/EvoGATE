"""
Formal Phase 2A human threshold sweep reproduction.

This runner supports:
- strict old-table targets from legacy GAT_results-string.csv
- additional requested feature combinations absent from the old table
- explicit no-degree Phase 2A runs
- consolidated TSV summaries under outputs/epgat_phase2a/human
"""

import argparse
import copy
from concurrent.futures import ProcessPoolExecutor, as_completed
import os

import pandas as pd
import yaml

from src.data.build_epgat_legacy_dataset import build_dataset
from src.train.train_epgat_legacy import train_from_config


def parse_args():
    parser = argparse.ArgumentParser(description="Run formal Phase 2A human threshold sweep")
    parser.add_argument("--config", required=True, type=str)
    return parser.parse_args()


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def write_yaml(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def _bool_cell(value):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in ["true", "1", "yes"]


def _combo_flags(combo_name):
    tokens = set([t for t in str(combo_name).split("_") if t])
    return {
        "expression": "EXP" in tokens,
        "orthologs": "ORT" in tokens,
        "sublocalization": "SUB" in tokens,
    }


def _load_old_table(config):
    df = pd.read_csv(config["target"]["source_file"])
    df = df[df["organism"].astype(str) == config["target"]["species"]].copy()
    df["feature_combo"] = (
        df["name"]
        .astype(str)
        .str.replace(config["target"].get("method_prefix", "GAT") + "_", "", n=1, regex=False)
    )
    return df


def _build_target_manifest(config, old_df):
    requested_combos = config["target"]["feature_combinations"]
    requested_thresholds = [int(v) for v in config["target"]["string_thresholds"]]
    combo_rank = {name: idx for idx, name in enumerate(requested_combos)}
    threshold_rank = {value: idx for idx, value in enumerate(requested_thresholds)}
    rows = []
    for combo_name in requested_combos:
        flags = _combo_flags(combo_name)
        combo_old = old_df[old_df["feature_combo"].astype(str) == combo_name].copy()
        for string_thr in requested_thresholds:
            matched = combo_old[combo_old["string_thr"].astype(int) == int(string_thr)].copy()
            row = {
                "feature_combo": combo_name,
                "row_name": "{}_{}".format(config["target"].get("method_prefix", "GAT"), combo_name),
                "string_thr": int(string_thr),
                "expression": bool(flags["expression"]),
                "orthologs": bool(flags["orthologs"]),
                "sublocalization": bool(flags["sublocalization"]),
                "old_result_available": not matched.empty,
                "old_auroc": "",
                "old_auprc": "",
                "old_accuracy": "",
                "old_mcc": "",
                "old_precision": "",
            }
            if not matched.empty:
                old_row = matched.iloc[0]
                row["old_auroc"] = float(old_row["mean"])
                row["old_auprc"] = float(old_row["auc_pr"])
                row["old_accuracy"] = float(old_row["accuracy"])
                row["old_mcc"] = float(old_row["mcc"]) if "mcc" in old_row.index else ""
                row["old_precision"] = float(old_row["precision"]) if "precision" in old_row.index else ""
            rows.append(row)
    manifest = pd.DataFrame(rows)
    manifest["combo_rank"] = manifest["feature_combo"].map(combo_rank)
    manifest["threshold_rank"] = manifest["string_thr"].astype(int).map(threshold_rank)
    manifest = manifest.sort_values(["combo_rank", "threshold_rank"]).reset_index(drop=True)
    return manifest.drop(columns=["combo_rank", "threshold_rank"])


def _aggregate_run_metrics(rows):
    frame = pd.DataFrame(rows)
    out = {}
    for metric in ["auroc", "auprc", "accuracy", "f1", "mcc", "best_val_auc", "test_count"]:
        out[metric + "_mean"] = float(frame[metric].mean())
        out[metric + "_std"] = float(frame[metric].std(ddof=0))
    return out


def _run_target_task(task_payload):
    config = task_payload["config"]
    target = task_payload["target"]
    combo_name = str(target["feature_combo"])
    string_thr = int(target["string_threshold"])
    task_run_root = os.path.join(config["paths"]["output_root"], config["run"]["name"], combo_name, "thr_{}".format(string_thr))
    os.makedirs(task_run_root, exist_ok=True)

    run_metrics = []
    for run_idx in range(int(config["execution"]["n_runs"])):
        run_seed = run_idx if config["execution"]["seed_strategy"] == "legacy_run_index" else int(config["base_config"]["legacy"]["seed"])
        run_config = copy.deepcopy(config["base_config"])
        run_config["run"]["name"] = os.path.join(config["run"]["name"], combo_name, "thr_{}".format(string_thr), "runs", "run_{}".format(run_idx))
        run_config["legacy"]["organism"] = config["target"]["species"]
        run_config["legacy"]["expression"] = bool(target["expression"])
        run_config["legacy"]["orthologs"] = bool(target["orthologs"])
        run_config["legacy"]["sublocalization"] = bool(target["sublocalization"])
        run_config["legacy"]["string_threshold"] = string_thr
        run_config["legacy"]["seed"] = run_seed
        build_dataset(run_config)
        result = train_from_config(run_config)
        row = dict(result["metrics"])
        row["run_id"] = "run_{}".format(run_idx)
        row["seed_used"] = run_seed
        run_metrics.append(row)

    all_runs = pd.DataFrame(run_metrics)
    all_runs.to_csv(os.path.join(task_run_root, "all_runs.tsv"), sep="\t", index=False)
    aggregated = _aggregate_run_metrics(run_metrics)
    task_result = {
        "feature_combo": combo_name,
        "row_name": str(target["row_name"]),
        "string_threshold": string_thr,
        "old_result_available": bool(target["old_result_available"]),
        "old_auroc": _to_float_or_blank(target["old_auroc"]),
        "old_auprc": _to_float_or_blank(target["old_auprc"]),
        "old_accuracy": _to_float_or_blank(target["old_accuracy"]),
        "old_mcc": _to_float_or_blank(target["old_mcc"]),
        "old_precision": _to_float_or_blank(target["old_precision"]),
        "auroc": float(aggregated["auroc_mean"]),
        "auprc": float(aggregated["auprc_mean"]),
        "mcc": float(aggregated["mcc_mean"]),
        "f1": float(aggregated["f1_mean"]),
        "acc": float(aggregated["accuracy_mean"]),
        "auroc_std": float(aggregated["auroc_std"]),
        "auprc_std": float(aggregated["auprc_std"]),
        "mcc_std": float(aggregated["mcc_std"]),
        "f1_std": float(aggregated["f1_std"]),
        "acc_std": float(aggregated["accuracy_std"]),
        "delta_auroc": "" if target["old_auroc"] == "" else float(aggregated["auroc_mean"]) - float(target["old_auroc"]),
        "delta_auprc": "" if target["old_auprc"] == "" else float(aggregated["auprc_mean"]) - float(target["old_auprc"]),
        "delta_acc": "" if target["old_accuracy"] == "" else float(aggregated["accuracy_mean"]) - float(target["old_accuracy"]),
        "n_runs": int(config["execution"]["n_runs"]),
        "epochs": int(config["base_config"]["train"]["epochs"]),
        "include_degree": bool(config["base_config"]["legacy"].get("include_degree", False)),
    }
    pd.DataFrame([task_result]).to_csv(os.path.join(task_run_root, "aggregated_metrics.tsv"), sep="\t", index=False)
    write_yaml(os.path.join(task_run_root, "run_config_frozen.yaml"), config)
    return task_result


def _to_float_or_blank(value):
    if value == "":
        return ""
    return float(value)


def main():
    args = parse_args()
    config = load_yaml(args.config)

    output_root = config["paths"]["output_root"]
    run_root = os.path.join(output_root, config["run"]["name"])
    os.makedirs(output_root, exist_ok=True)
    os.makedirs(run_root, exist_ok=True)

    old_df = _load_old_table(config)
    target_manifest = _build_target_manifest(config, old_df)
    target_manifest.to_csv(os.path.join(output_root, "target_manifest.tsv"), sep="\t", index=False)

    threshold_rows = []
    pending_payloads = []
    for _, target in target_manifest.iterrows():
        combo_name = str(target["feature_combo"])
        string_thr = int(target["string_thr"])
        target_run_root = os.path.join(run_root, combo_name, "thr_{}".format(string_thr))
        aggregated_path = os.path.join(target_run_root, "aggregated_metrics.tsv")
        if bool(config.get("execution", {}).get("resume_existing", True)) and os.path.exists(aggregated_path):
            existing = pd.read_csv(aggregated_path, sep="\t").iloc[0].to_dict()
            threshold_rows.append(existing)
            continue
        target_payload = target.to_dict()
        target_payload["string_threshold"] = string_thr
        pending_payloads.append({"config": config, "target": target_payload})

    max_parallel = int(config.get("execution", {}).get("max_parallel", 1))
    if max_parallel <= 1:
        for payload in pending_payloads:
            threshold_rows.append(_run_target_task(payload))
    else:
        with ProcessPoolExecutor(max_workers=max_parallel) as pool:
            futures = [pool.submit(_run_target_task, payload) for payload in pending_payloads]
            for future in as_completed(futures):
                threshold_rows.append(future.result())

    threshold_summary = pd.DataFrame(threshold_rows).sort_values(["feature_combo", "string_threshold"]).reset_index(drop=True)
    threshold_summary.to_csv(os.path.join(output_root, "threshold_sweep_summary.tsv"), sep="\t", index=False)

    feature_rows = []
    for combo_name, combo_df in threshold_summary.groupby("feature_combo", sort=True):
        best_idx = combo_df["auroc"].astype(float).idxmax()
        best_row = threshold_summary.loc[best_idx]
        old_df_combo = combo_df[combo_df["old_result_available"] == True].copy()
        old_best_idx = old_df_combo["old_auroc"].astype(float).idxmax() if not old_df_combo.empty else None
        old_best_auroc = threshold_summary.loc[old_best_idx, "old_auroc"] if old_best_idx is not None else ""
        feature_rows.append(
            {
                "feature_combo": combo_name,
                "best_threshold": int(best_row["string_threshold"]),
                "best_auroc": float(best_row["auroc"]),
                "best_auprc": float(best_row["auprc"]),
                "best_mcc": float(best_row["mcc"]),
                "best_f1": float(best_row["f1"]),
                "best_acc": float(best_row["acc"]),
                "old_result_available": bool(combo_df["old_result_available"].any()),
                "best_old_auroc": old_best_auroc,
                "delta_vs_best_old_auroc": "" if old_best_auroc == "" else float(best_row["auroc"]) - float(old_best_auroc),
            }
        )
    feature_summary = pd.DataFrame(feature_rows).sort_values(["best_auroc", "feature_combo"], ascending=[False, True]).reset_index(drop=True)
    feature_summary.to_csv(os.path.join(output_root, "feature_combo_summary.tsv"), sep="\t", index=False)

    old_available = threshold_summary[threshold_summary["old_result_available"] == True].copy()
    best_old_match_idx = old_available["delta_auroc"].abs().astype(float).idxmin() if not old_available.empty else None
    best_reproduced_idx = threshold_summary["auroc"].astype(float).idxmax()
    best_reproduced = threshold_summary.loc[best_reproduced_idx]
    best_old_match = threshold_summary.loc[best_old_match_idx] if best_old_match_idx is not None else None

    final_rows = [
        {
            "summary_type": "best_reproduced_overall",
            "feature_combo": best_reproduced["feature_combo"],
            "string_threshold": int(best_reproduced["string_threshold"]),
            "auroc": float(best_reproduced["auroc"]),
            "auprc": float(best_reproduced["auprc"]),
            "mcc": float(best_reproduced["mcc"]),
            "f1": float(best_reproduced["f1"]),
            "acc": float(best_reproduced["acc"]),
            "old_result_available": bool(best_reproduced["old_result_available"]),
            "old_auroc": _to_float_or_blank(best_reproduced["old_auroc"]),
            "delta_auroc": _to_float_or_blank(best_reproduced["delta_auroc"]),
        }
    ]
    if best_old_match is not None:
        final_rows.append(
            {
                "summary_type": "closest_old_match",
                "feature_combo": best_old_match["feature_combo"],
                "string_threshold": int(best_old_match["string_threshold"]),
                "auroc": float(best_old_match["auroc"]),
                "auprc": float(best_old_match["auprc"]),
                "mcc": float(best_old_match["mcc"]),
                "f1": float(best_old_match["f1"]),
                "acc": float(best_old_match["acc"]),
                "old_result_available": bool(best_old_match["old_result_available"]),
                "old_auroc": float(best_old_match["old_auroc"]),
                "delta_auroc": float(best_old_match["delta_auroc"]),
            }
        )
    final_summary = pd.DataFrame(final_rows)
    final_summary.to_csv(os.path.join(output_root, "final_reproduced_summary.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
