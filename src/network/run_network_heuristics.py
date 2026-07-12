import argparse
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from src.classical_baselines.common import (
    HEURISTIC_METHODS,
    build_dataset_for_benchmark,
    build_prediction_table,
    compute_binary_metrics,
    dump_yaml,
    load_dataset_bundle,
    load_yaml,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Run deterministic network heuristics for the classical baseline benchmark")
    parser.add_argument("--config", required=True, type=str)
    parser.add_argument("--species", required=True, type=str)
    parser.add_argument("--method", required=True, choices=sorted(HEURISTIC_METHODS))
    parser.add_argument("--output-dir", required=True, type=str)
    return parser.parse_args()


def compute_scores(method, edge_index, num_nodes):
    graph = nx.Graph()
    graph.add_nodes_from(range(int(num_nodes)))
    graph.add_edges_from((int(src), int(dst)) for src, dst in edge_index.tolist())
    if method == "DC":
        scores_map = nx.degree_centrality(graph)
        score_name = "degree_centrality"
    elif method == "CC":
        scores_map = nx.clustering(graph)
        score_name = "clustering_coefficient"
    else:
        raise ValueError(f"Unsupported heuristic method: {method}")
    scores = np.array([float(scores_map.get(node_idx, 0.0)) for node_idx in range(int(num_nodes))], dtype=np.float64)
    return scores, score_name


def top_k_predictions_for_labeled(scores, labeled_indices, positive_total):
    pred = np.zeros(scores.shape[0], dtype=np.int64)
    if positive_total <= 0 or labeled_indices.size == 0:
        return pred
    ranked = labeled_indices[np.argsort(-scores[labeled_indices], kind="mergesort")]
    pred[ranked[: int(positive_total)]] = 1
    return pred


def main():
    args = parse_args()
    config = load_yaml(args.config)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset_meta = build_dataset_for_benchmark(
        config=config,
        species=args.species,
        feature_setting="network",
        output_dir=output_dir,
        seed=int(config["runtime"]["base_seed"]),
    )
    bundle = load_dataset_bundle(output_dir)

    scores, score_name = compute_scores(args.method, bundle["edge_index"], len(bundle["node_manifest"]))
    labeled_indices = np.array(
        [bundle["mapping"][gene_id] for gene_id in bundle["labeled"]["legacy_gene_id"].astype(str)],
        dtype=np.int64,
    )
    positive_total = int((bundle["labeled"]["label_numeric"].astype(float) == 1).sum())
    pred_all = top_k_predictions_for_labeled(scores, labeled_indices, positive_total)

    test_idx = bundle["test_idx"]
    y_test = bundle["y_all"][test_idx].astype(int)
    test_scores = scores[test_idx]
    test_pred = pred_all[test_idx].astype(int)
    metrics = compute_binary_metrics(y_test, test_scores, test_pred)

    metrics_row = {
        "species": args.species,
        "method": args.method,
        "feature_setting": "network",
        "label_regime": dataset_meta["label_regime"],
        "score_name": score_name,
        "run_id": "network",
        "test_count": int(len(test_idx)),
        **metrics,
    }
    pd.DataFrame([metrics_row]).to_csv(output_dir / "metrics.tsv", sep="\t", index=False)

    prediction_table = build_prediction_table(
        bundle["node_manifest"],
        bundle["label_manifest"],
        pred_score=scores,
        pred_label=pred_all,
        method=args.method,
        feature_setting="network",
    )
    prediction_table.to_csv(output_dir / "predictions.tsv", sep="\t", index=False)

    method_summary = {
        "species": args.species,
        "method": args.method,
        "feature_setting": "network",
        "label_regime": dataset_meta["label_regime"],
        "node_count": int(len(bundle["node_manifest"])),
        "edge_count": int(bundle["edge_index"].shape[0]),
        "labeled_count": int(len(bundle["labeled"])),
        "positive_labeled_count": int(positive_total),
        "test_count": int(len(test_idx)),
        "score_name": score_name,
        "prediction_rule": "top_k_equals_total_positive_labeled",
        "positive_set_path": dataset_meta["positive_set_path"],
        "negative_set_path": dataset_meta["negative_set_path"],
    }
    pd.DataFrame([method_summary]).to_csv(output_dir / "method_summary.tsv", sep="\t", index=False)

    resolved_config = {
        "species": args.species,
        "method": args.method,
        "feature_setting": "network",
        "output_dir": str(output_dir),
        "label_regime": dataset_meta["label_regime"],
        "dataset_config": dataset_meta["builder_config"],
        "score_name": score_name,
        "prediction_rule": "top_k_equals_total_positive_labeled",
        "metrics": metrics_row,
    }
    dump_yaml(resolved_config, output_dir / "resolved_config.yaml")


if __name__ == "__main__":
    main()
