import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import yaml

from src.data.freeze_unified_protocol import load_config
from src.eval.aggregate_frozen_protocol_runs import aggregate_runs
from src.train.run_frozen_protocol_model import DETERMINISTIC_MODELS, run_benchmark_task


DEFAULT_PROTOCOLS = ["fgraminearum_oldlabel", "fgraminearum_newlabel", "scerevisiae"]
DEFAULT_MODELS = ["N2V_MLP", "GAT", "GCN", "GIN", "GraphSAGE", "MLP", "DC", "CC"]


def parse_args():
    parser = argparse.ArgumentParser(description="Run the clean true-Node2Vec frozen benchmark subset")
    parser.add_argument("--config", default="configs/frozen_protocol.yaml", type=str)
    parser.add_argument(
        "--output-root",
        default="results/frozen_protocol_benchmark_v3_true_n2v",
        type=str,
    )
    parser.add_argument("--protocols", nargs="*", default=DEFAULT_PROTOCOLS)
    parser.add_argument("--models", nargs="*", default=DEFAULT_MODELS)
    return parser.parse_args()


def _assert_true_n2v_runtime(config):
    if str(config["runtime"].get("node2vec_backend", "")).strip().lower() == "svd":
        raise RuntimeError("Clean subset benchmark refuses node2vec_backend=svd")
    if not bool(config["runtime"].get("require_true_node2vec", False)):
        raise RuntimeError("Clean subset benchmark requires runtime.require_true_node2vec=true")


def _write_clean_subset_config(source_config_path, output_root):
    config = load_config(source_config_path)
    config["runtime"]["seed_list"] = [1029]
    config["runtime"]["node2vec_backend"] = "native_walk"
    config["runtime"]["require_true_node2vec"] = True
    config["runtime"]["graph_contract"] = "undirected_symmetrized"
    config["models"]["N2V_MLP"]["walk_length"] = 8
    config["models"]["N2V_MLP"]["context_size"] = 4
    config["models"]["N2V_MLP"]["walks_per_node"] = 2
    config["models"]["N2V_MLP"]["epochs"] = 2
    config["models"]["N2V_MLP"]["batch_size"] = 1024
    config["models"]["N2V_MLP"]["p"] = 1.0
    config["models"]["N2V_MLP"]["q"] = 1.0
    clean_config_path = output_root / "clean_subset_config.yaml"
    clean_config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return clean_config_path, config


def _resolve_clean_subset_config(source_config_path, output_root):
    clean_config_path = output_root / "clean_subset_config.yaml"
    if clean_config_path.exists():
        config = load_config(clean_config_path)
        _assert_true_n2v_runtime(config)
        if str(config["runtime"].get("graph_contract", "")).strip().lower() != "undirected_symmetrized":
            raise RuntimeError("Existing clean subset config must keep graph_contract=undirected_symmetrized")
        return clean_config_path, config
    return _write_clean_subset_config(source_config_path, output_root)


def _task_output_dir(config, output_root, protocol, model, seed):
    feature_setting = str(config["models"][model]["feature_setting"])
    if model in DETERMINISTIC_MODELS:
        return output_root / "runs" / protocol / model / feature_setting / "deterministic"
    return output_root / "runs" / protocol / model / feature_setting / f"run_{seed}"


def _collect_existing_run(task_output):
    metrics_path = task_output / "metrics.tsv"
    if not metrics_path.exists():
        return None
    row = pd.read_csv(metrics_path, sep="\t").iloc[0].to_dict()
    row["_status"] = "completed_existing"
    row["_output_dir"] = str(task_output)
    return row


def _run_tasks(config_path, output_root, protocols, models):
    config = load_config(config_path)
    _assert_true_n2v_runtime(config)
    seeds = [int(seed) for seed in config["runtime"]["seed_list"]]
    rows = []
    status_rows = []
    run_root = output_root / "runs"
    run_root.mkdir(parents=True, exist_ok=True)
    for protocol in protocols:
        for model in models:
            if model in DETERMINISTIC_MODELS:
                task_output = _task_output_dir(config, output_root, protocol, model, None)
                existing = _collect_existing_run(task_output)
                if existing is not None:
                    rows.append(existing)
                    status_rows.append(
                        {
                            "protocol": protocol,
                            "model": model,
                            "feature_setting": str(config["models"][model]["feature_setting"]),
                            "run_id": "deterministic",
                            "seed": "",
                            "expected_output_dir": str(task_output),
                            "metrics_present": "true",
                            "status": "completed_existing",
                        }
                    )
                    continue
                result = run_benchmark_task(
                    config_path=config_path,
                    protocol_name=protocol,
                    model_name=model,
                    output_dir=task_output,
                    seed=None,
                    graph_contract=config["runtime"]["graph_contract"],
                )
                result["_status"] = "completed_new"
                result["_output_dir"] = str(task_output)
                rows.append(result)
                status_rows.append(
                    {
                        "protocol": protocol,
                        "model": model,
                        "feature_setting": str(config["models"][model]["feature_setting"]),
                        "run_id": "deterministic",
                        "seed": "",
                        "expected_output_dir": str(task_output),
                        "metrics_present": "true",
                        "status": "completed_new",
                    }
                )
            else:
                for seed in seeds:
                    task_output = _task_output_dir(config, output_root, protocol, model, seed)
                    existing = _collect_existing_run(task_output)
                    if existing is not None:
                        rows.append(existing)
                        status_rows.append(
                            {
                                "protocol": protocol,
                                "model": model,
                                "feature_setting": str(config["models"][model]["feature_setting"]),
                                "run_id": f"seed_{seed}",
                                "seed": int(seed),
                                "expected_output_dir": str(task_output),
                                "metrics_present": "true",
                                "status": "completed_existing",
                            }
                        )
                        continue
                    result = run_benchmark_task(
                        config_path=config_path,
                        protocol_name=protocol,
                        model_name=model,
                        output_dir=task_output,
                        seed=seed,
                        graph_contract=config["runtime"]["graph_contract"],
                    )
                    result["_status"] = "completed_new"
                    result["_output_dir"] = str(task_output)
                    rows.append(result)
                    status_rows.append(
                        {
                            "protocol": protocol,
                            "model": model,
                            "feature_setting": str(config["models"][model]["feature_setting"]),
                            "run_id": f"seed_{seed}",
                            "seed": int(seed),
                            "expected_output_dir": str(task_output),
                            "metrics_present": "true",
                            "status": "completed_new",
                        }
                    )
    return pd.DataFrame(rows), pd.DataFrame(status_rows)


def _build_split_manifest(per_run):
    rows = []
    for protocol, group in per_run.groupby("protocol", sort=True):
        split_path = Path(str(group["split_manifest"].iloc[0]))
        split_df = pd.read_csv(split_path, sep="\t")
        split_df["label"] = pd.to_numeric(split_df["label"], errors="raise").astype(int)
        for split_name in ["train", "val", "test"]:
            mask = split_df["split"].astype(str) == split_name
            rows.append(
                {
                    "protocol": protocol,
                    "split_version": str(group["split_version"].iloc[0]),
                    "split_name": split_name,
                    "count": int(mask.sum()),
                    "positives": int((mask & (split_df["label"] == 1)).sum()),
                    "negatives": int((mask & (split_df["label"] == 0)).sum()),
                }
            )
    return pd.DataFrame(rows).sort_values(["protocol", "split_name"]).reset_index(drop=True)


def _build_label_regime_manifest(config, protocols):
    rows = []
    for protocol in protocols:
        protocol_cfg = config["protocols"][protocol]
        label_path = Path(config["paths"]["labels_dir"]) / str(protocol_cfg["label_output"])
        label_df = pd.read_csv(label_path, sep="\t")
        rows.append(
            {
                "protocol": protocol,
                "species": str(protocol_cfg["species"]),
                "label_regime": str(protocol_cfg["regime"]),
                "label_manifest": str(label_path),
                "label_count": int(len(label_df)),
                "positive_count": int((pd.to_numeric(label_df["label"], errors="raise").astype(int) == 1).sum()),
                "negative_count": int((pd.to_numeric(label_df["label"], errors="raise").astype(int) == 0).sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(["species", "label_regime"]).reset_index(drop=True)


def _build_graph_contract_manifest(per_run):
    columns = [
        "protocol",
        "species",
        "label_regime",
        "model",
        "feature_setting",
        "feature_contract_group",
        "graph_contract",
        "graph_source",
        "edge_count",
        "node_count",
        "split_version",
        "threshold_strategy",
        "evaluation_contract",
    ]
    return per_run[columns].drop_duplicates().sort_values(["protocol", "model"]).reset_index(drop=True)


def _build_node2vec_backend_manifest(output_root):
    rows = []
    for summary_path in sorted(output_root.glob("runs/*/N2V_MLP/N2V/run_*/node2vec_summary.tsv")):
        row = pd.read_csv(summary_path, sep="\t").iloc[0].to_dict()
        row["summary_path"] = str(summary_path)
        rows.append(row)
    return pd.DataFrame(rows)


def _write_run_status(run_status, output_path):
    ordered = [
        "protocol",
        "model",
        "feature_setting",
        "run_id",
        "seed",
        "expected_output_dir",
        "metrics_present",
        "status",
    ]
    if run_status.empty:
        pd.DataFrame(columns=ordered).to_csv(output_path, sep="\t", index=False)
        return
    run_status = run_status[ordered].sort_values(["protocol", "model", "run_id"], kind="stable").reset_index(drop=True)
    run_status.to_csv(output_path, sep="\t", index=False)


def _build_sage_aggregator_check():
    try:
        from src.models.epgat_sage import EPGATOriginalSAGE
    except ImportError as exc:
        return pd.DataFrame(
            [
                {
                    "aggregator_left": "mean",
                    "aggregator_right": "pool",
                    "same_shape": "na",
                    "allclose": "na",
                    "l2_distance": np.nan,
                    "max_abs_diff": np.nan,
                    "status": f"skipped_import_error: {exc}",
                }
            ]
        )
    torch.manual_seed(1029)
    x = torch.randn(5, 4)
    edge_index = torch.tensor(
        [
            [0, 0, 1, 2, 3, 4, 4],
            [1, 2, 2, 3, 4, 0, 1],
        ],
        dtype=torch.long,
    )
    mean_model = EPGATOriginalSAGE(in_feats=4, n_hidden=8, n_layers=2, dropout=0.0, aggregator_type="mean")
    torch.manual_seed(1029)
    pool_model = EPGATOriginalSAGE(in_feats=4, n_hidden=8, n_layers=2, dropout=0.0, aggregator_type="pool")
    mean_model.eval()
    pool_model.eval()
    with torch.no_grad():
        mean_out = mean_model(x, edge_index)
        pool_out = pool_model(x, edge_index)
    diff = pool_out - mean_out
    return pd.DataFrame(
        [
            {
                "aggregator_left": "mean",
                "aggregator_right": "pool",
                "same_shape": str(tuple(mean_out.shape) == tuple(pool_out.shape)).lower(),
                "allclose": str(bool(torch.allclose(mean_out, pool_out))).lower(),
                "l2_distance": float(torch.norm(diff).item()),
                "max_abs_diff": float(diff.abs().max().item()),
                "status": "ok",
            }
        ]
    )


def _summarize_benchmark(aggregated):
    summary = aggregated.copy()
    summary["leaderboard_group"] = summary.apply(
        lambda row: "centrality_top_k"
        if str(row["threshold_strategy"]) == "top_k_labeled_positive_count"
        else str(row["feature_contract_group"]),
        axis=1,
    )
    summary["rank_within_group"] = (
        summary.groupby(["protocol", "leaderboard_group"], sort=True)["test_auprc_mean"]
        .rank(ascending=False, method="dense")
        .astype(int)
    )
    keep_columns = [
        "protocol",
        "species",
        "label_regime",
        "model",
        "feature_setting",
        "leaderboard_group",
        "graph_contract",
        "split_version",
        "threshold_strategy",
        "evaluation_contract",
        "embedding_backend",
        "fallback_used",
        "n_runs",
        "test_auprc_mean",
        "test_auroc_mean",
        "test_mcc_mean",
        "test_specificity_mean",
        "rank_within_group",
    ]
    return summary[keep_columns].sort_values(
        ["protocol", "leaderboard_group", "rank_within_group", "test_auprc_mean"],
        ascending=[True, True, True, False],
    ).reset_index(drop=True)


def _write_markdown(summary_df, output_path):
    columns = list(summary_df.columns)
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = [header, divider]
    for _, row in summary_df.iterrows():
        rows.append("| " + " | ".join(str(row[column]) for column in columns) + " |")
    lines = [
        "# True Node2Vec Clean Benchmark Summary",
        "",
        "Separate leaderboards are reported for `topology_embedding_contract`, `ort_exp_sub_contract`, and `centrality_top_k`.",
        "",
        *rows,
        "",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")


def _write_fairness_assertions(per_run, output_path):
    fair_models = per_run[~per_run["model"].isin(["DC", "CC"])].copy()
    graph_contracts = sorted(set(fair_models["graph_contract"].astype(str)))
    thresholds = sorted(set(fair_models["threshold_strategy"].astype(str)))
    eval_contracts = sorted(set(fair_models["evaluation_contract"].astype(str)))
    n2v_rows = fair_models[fair_models["model"] == "N2V_MLP"].copy()
    fake_rows = n2v_rows[n2v_rows["fallback_used"].astype(str).str.lower() == "true"]
    lines = [
        "# Fairness Assertions",
        "",
        f"- Shared graph contract across trainable fair-comparison models: `{', '.join(graph_contracts)}`.",
        f"- Shared threshold strategy across trainable fair-comparison models: `{', '.join(thresholds)}`.",
        f"- Shared evaluation contract across trainable fair-comparison models: `{', '.join(eval_contracts)}`.",
        "- Fusarium `oldlabel` and `newlabel` remain separate protocols and are not merged into one leaderboard.",
        "- `DC` and `CC` keep `top_k_labeled_positive_count` and are therefore reported in a separate `centrality_top_k` leaderboard.",
        f"- Node2Vec fallback rows detected in this clean subset: `{int(len(fake_rows))}`.",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    clean_config_path, config = _resolve_clean_subset_config(args.config, output_root)
    per_run, run_status = _run_tasks(str(clean_config_path), output_root, args.protocols, args.models)
    per_run = per_run.drop(columns=[column for column in ["_status", "_output_dir"] if column in per_run.columns])
    aggregated = aggregate_runs(per_run)
    benchmark_summary = _summarize_benchmark(aggregated)
    graph_contract_manifest = _build_graph_contract_manifest(per_run)
    split_manifest = _build_split_manifest(per_run)
    label_regime_manifest = _build_label_regime_manifest(config, args.protocols)
    node2vec_backend_manifest = _build_node2vec_backend_manifest(output_root)
    sage_aggregator_check = _build_sage_aggregator_check()

    benchmark_summary.to_csv(output_root / "benchmark_summary_true_n2v.tsv", sep="\t", index=False)
    _write_markdown(benchmark_summary, output_root / "benchmark_summary_true_n2v.md")
    graph_contract_manifest.to_csv(output_root / "graph_contract_manifest.tsv", sep="\t", index=False)
    split_manifest.to_csv(output_root / "split_manifest.tsv", sep="\t", index=False)
    label_regime_manifest.to_csv(output_root / "label_regime_manifest.tsv", sep="\t", index=False)
    node2vec_backend_manifest.to_csv(output_root / "node2vec_backend_manifest.tsv", sep="\t", index=False)
    sage_aggregator_check.to_csv(output_root / "sage_aggregator_check.tsv", sep="\t", index=False)
    _write_fairness_assertions(per_run, output_root / "fairness_assertions.md")
    _write_run_status(run_status, output_root / "run_status.tsv")


if __name__ == "__main__":
    main()
