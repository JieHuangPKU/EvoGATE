import argparse
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


METRIC_COLS = ["test_auprc", "test_mcc", "val_auprc", "val_mcc", "test_auroc", "test_accuracy"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize Figure 4 Graph Robustness outputs")
    parser.add_argument("--protocol", required=True, type=str)
    parser.add_argument("--model", required=True, type=str)
    parser.add_argument("--feature-setting", required=True, type=str)
    parser.add_argument("--thresholds", required=True, type=str)
    parser.add_argument("--seeds", required=True, type=str)
    parser.add_argument("--threshold-output-root", required=True, type=str)
    parser.add_argument("--main-efg-output-root", required=True, type=str)
    parser.add_argument("--supplementary-output-roots", default="", type=str)
    parser.add_argument("--string-graph", required=True, type=str)
    parser.add_argument("--main-efg-graph", required=True, type=str)
    parser.add_argument("--supplementary-efg-graphs", default="", type=str)
    parser.add_argument("--split-manifest", required=True, type=str)
    parser.add_argument("--efg-adapter-summary", required=True, type=str)
    parser.add_argument("--supplementary-efg-adapter-summaries", default="", type=str)
    parser.add_argument("--summary-dir", required=True, type=str)
    return parser.parse_args()


def parse_int_csv(text: str) -> list[int]:
    return [int(token.strip()) for token in str(text).split(",") if token.strip()]


def parse_mapping_csv(text: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if not str(text).strip():
        return mapping
    for token in str(text).split(","):
        token = token.strip()
        if not token:
            continue
        if "=" not in token:
            raise ValueError(f"Expected NAME=VALUE entry, got '{token}'")
        key, value = token.split("=", 1)
        mapping[key.strip()] = value.strip()
    return mapping


def load_graph_edges(path: Path, threshold: int | None = None) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python", dtype=str).fillna("")
    if "combined_score" in df.columns and threshold is not None:
        score = pd.to_numeric(df["combined_score"], errors="coerce").fillna(0.0)
        df = df.loc[score >= float(threshold)].copy()
    if "A" not in df.columns or "B" not in df.columns:
        raise ValueError(f"Graph file must contain A/B columns: {path}")
    df = df[["A", "B"]].copy()
    df["A"] = df["A"].astype(str).str.strip()
    df["B"] = df["B"].astype(str).str.strip()
    df = df[df["A"].ne("") & df["B"].ne("") & df["A"].ne(df["B"])].copy()
    df[["node_u", "node_v"]] = pd.DataFrame(
        np.sort(df[["A", "B"]].to_numpy(dtype=object), axis=1),
        index=df.index,
    )
    return df.drop_duplicates(["node_u", "node_v"]).reset_index(drop=True)


def undirected_edge_set(edge_df: pd.DataFrame) -> set[tuple[str, str]]:
    return set(edge_df[["node_u", "node_v"]].itertuples(index=False, name=None))


def compute_graph_stats(edge_df: pd.DataFrame, labeled_nodes: set[str], graph_label: str, threshold: int | None, model: str, feature_setting: str) -> dict[str, object]:
    graph_nodes = set(edge_df["node_u"].astype(str)) | set(edge_df["node_v"].astype(str))
    degree_map = {node: 0 for node in labeled_nodes}
    for node_u, node_v in edge_df[["node_u", "node_v"]].itertuples(index=False, name=None):
        if node_u in degree_map:
            degree_map[node_u] += 1
        if node_v in degree_map:
            degree_map[node_v] += 1

    graph = nx.Graph()
    graph.add_nodes_from(graph_nodes)
    graph.add_edges_from(edge_df[["node_u", "node_v"]].itertuples(index=False, name=None))
    largest_cc = max((len(component) for component in nx.connected_components(graph)), default=0)

    return {
        "graph_condition": graph_label,
        "graph_family": "STRING",
        "string_threshold": threshold if threshold is not None else "",
        "model": model,
        "feature_setting": feature_setting,
        "model_feature_setting": f"{model} + {feature_setting}",
        "edge_count": int(len(edge_df)),
        "graph_node_count": int(len(graph_nodes)),
        "labeled_node_count": int(len(labeled_nodes)),
        "isolated_node_count": int(sum(1 for value in degree_map.values() if value == 0)),
        "mean_degree_labeled_nodes": float(np.mean(list(degree_map.values())) if degree_map else 0.0),
        "largest_connected_component_size": int(largest_cc),
    }


def load_run_metrics(path: Path) -> dict[str, object]:
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        raise ValueError(f"Metrics file is empty: {path}")
    return df.iloc[0].to_dict()


def aggregate_metrics(df: pd.DataFrame, group_cols: list[str], metric_cols: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for keys, group_df in df.groupby(group_cols, dropna=False, sort=False):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        row = {col: value for col, value in zip(group_cols, key_values)}
        row["n_runs"] = int(len(group_df))
        for metric in metric_cols:
            values = pd.to_numeric(group_df[metric], errors="coerce")
            row[f"{metric}_mean"] = float(values.mean())
            row[f"{metric}_sd"] = float(values.std(ddof=0))
            row[f"{metric}_sem"] = float(values.std(ddof=0) / np.sqrt(len(values))) if len(values) > 0 else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def ordered_source_summary(source_df: pd.DataFrame, order: list[str]) -> pd.DataFrame:
    if source_df.empty:
        return source_df
    return source_df.set_index("graph_condition").loc[order].reset_index()


def annotate_metrics(metrics: dict[str, object], graph_condition: str, graph_family: str, string_threshold: object, run_group: str, metrics_path: Path, model: str, feature_setting: str) -> dict[str, object]:
    metrics = dict(metrics)
    metrics["graph_condition"] = graph_condition
    metrics["graph_family"] = graph_family
    metrics["string_threshold"] = string_threshold
    metrics["run_group"] = run_group
    metrics["metrics_path"] = str(metrics_path)
    metrics.setdefault("model", model)
    metrics.setdefault("feature_setting", feature_setting)
    metrics["model_feature_setting"] = f"{model} + {feature_setting}"
    return metrics


def main() -> None:
    args = parse_args()
    thresholds = parse_int_csv(args.thresholds)
    seeds = parse_int_csv(args.seeds)
    threshold_output_root = Path(args.threshold_output_root)
    main_efg_output_root = Path(args.main_efg_output_root)
    supplementary_output_roots = parse_mapping_csv(args.supplementary_output_roots)
    string_graph = Path(args.string_graph)
    main_efg_graph = Path(args.main_efg_graph)
    supplementary_efg_graphs = parse_mapping_csv(args.supplementary_efg_graphs)
    supplementary_adapter_summaries = parse_mapping_csv(args.supplementary_efg_adapter_summaries)
    split_manifest = pd.read_csv(args.split_manifest, sep="\t", dtype=str).fillna("")
    summary_dir = Path(args.summary_dir)
    summary_dir.mkdir(parents=True, exist_ok=True)

    labeled_nodes = set(split_manifest["graph_gene_id"].astype(str))

    per_run_rows: list[dict[str, object]] = []
    for threshold in thresholds:
        for seed in seeds:
            metrics_path = threshold_output_root / f"string_thr_{threshold}" / args.protocol / args.model / args.feature_setting / f"run_{seed}" / "metrics.tsv"
            metrics = load_run_metrics(metrics_path)
            per_run_rows.append(
                annotate_metrics(
                    metrics,
                    graph_condition=f"STRING_{threshold}",
                    graph_family="STRING",
                    string_threshold=threshold,
                    run_group="threshold_sweep",
                    metrics_path=metrics_path,
                    model=args.model,
                    feature_setting=args.feature_setting,
                )
            )

    for seed in seeds:
        metrics_path = main_efg_output_root / args.protocol / args.model / args.feature_setting / f"run_{seed}" / "metrics.tsv"
        metrics = load_run_metrics(metrics_path)
        per_run_rows.append(
            annotate_metrics(
                metrics,
                graph_condition="eFG",
                graph_family="eFG",
                string_threshold="",
                run_group="source_comparison",
                metrics_path=metrics_path,
                model=args.model,
                feature_setting=args.feature_setting,
            )
        )

    per_run_df = pd.DataFrame(per_run_rows)
    per_run_df.to_csv(summary_dir / "Figure4_threshold_per_run_metrics.tsv", sep="\t", index=False)

    threshold_df = per_run_df[per_run_df["graph_family"].astype(str).eq("STRING")].copy()
    threshold_agg_df = aggregate_metrics(
        threshold_df,
        ["graph_condition", "graph_family", "string_threshold", "model", "feature_setting", "model_feature_setting"],
        METRIC_COLS,
    )
    threshold_agg_df = threshold_agg_df.sort_values("string_threshold").reset_index(drop=True)
    threshold_agg_df.to_csv(summary_dir / "Figure4_threshold_aggregated_metrics.tsv", sep="\t", index=False)

    density_rows = [
        compute_graph_stats(load_graph_edges(string_graph, threshold=threshold), labeled_nodes, f"STRING_{threshold}", threshold, args.model, args.feature_setting)
        for threshold in thresholds
    ]
    density_df = pd.DataFrame(density_rows).sort_values("string_threshold").reset_index(drop=True)
    density_df.to_csv(summary_dir / "Figure4_network_density_summary.tsv", sep="\t", index=False)

    source_per_run_df = per_run_df[per_run_df["graph_condition"].isin(["STRING_300", "STRING_700", "eFG"])].copy()
    source_per_run_df.to_csv(summary_dir / "Figure4_source_comparison_per_run_metrics.tsv", sep="\t", index=False)
    source_metrics_df = aggregate_metrics(
        source_per_run_df,
        ["graph_condition", "graph_family", "model", "feature_setting", "model_feature_setting"],
        METRIC_COLS,
    )
    source_metrics_df["error_bar_definition"] = "mean_plus_minus_sd"
    source_metrics_df = ordered_source_summary(source_metrics_df, ["STRING_300", "STRING_700", "eFG"])
    source_metrics_df.to_csv(summary_dir / "Figure4_source_comparison_metrics.tsv", sep="\t", index=False)

    string_300_edges = load_graph_edges(string_graph, threshold=300)
    string_700_edges = load_graph_edges(string_graph, threshold=700)
    efg_edges = load_graph_edges(main_efg_graph)
    set_300 = undirected_edge_set(string_300_edges)
    set_700 = undirected_edge_set(string_700_edges)
    set_efg = undirected_edge_set(efg_edges)
    overlap_row = {
        "protocol": args.protocol,
        "model": args.model,
        "feature_setting": args.feature_setting,
        "model_feature_setting": f"{args.model} + {args.feature_setting}",
        "efg_main_version": "eFG_ALL",
        "string_300_edge_count": int(len(set_300)),
        "string_700_edge_count": int(len(set_700)),
        "efg_edge_count": int(len(set_efg)),
        "string_300_intersect_efg_edge_count": int(len(set_300 & set_efg)),
        "string_700_intersect_efg_edge_count": int(len(set_700 & set_efg)),
        "string_300_only_edge_count": int(len(set_300 - set_efg)),
        "string_700_only_edge_count": int(len(set_700 - set_efg)),
        "efg_only_vs_string_300_edge_count": int(len(set_efg - set_300)),
        "efg_only_vs_string_700_edge_count": int(len(set_efg - set_700)),
        "string_300_overlap_ratio": float(len(set_300 & set_efg) / len(set_300)) if set_300 else 0.0,
        "string_700_overlap_ratio": float(len(set_700 & set_efg) / len(set_700)) if set_700 else 0.0,
        "efg_overlap_ratio_vs_string_300": float(len(set_300 & set_efg) / len(set_efg)) if set_efg else 0.0,
        "efg_overlap_ratio_vs_string_700": float(len(set_700 & set_efg) / len(set_efg)) if set_efg else 0.0,
        "string_300_node_overlap_with_efg": int(len((set(string_300_edges["node_u"]) | set(string_300_edges["node_v"])) & (set(efg_edges["node_u"]) | set(efg_edges["node_v"])))),
        "string_700_node_overlap_with_efg": int(len((set(string_700_edges["node_u"]) | set(string_700_edges["node_v"])) & (set(efg_edges["node_u"]) | set(efg_edges["node_v"])))),
        "efg_adapter_summary_path": str(Path(args.efg_adapter_summary)),
    }
    pd.DataFrame([overlap_row]).to_csv(summary_dir / "Figure4_edge_overlap_summary.tsv", sep="\t", index=False)

    supplementary_rows: list[dict[str, object]] = []
    for graph_condition, output_root in supplementary_output_roots.items():
        for seed in seeds:
            metrics_path = Path(output_root) / args.protocol / args.model / args.feature_setting / f"run_{seed}" / "metrics.tsv"
            metrics = load_run_metrics(metrics_path)
            supplementary_rows.append(
                annotate_metrics(
                    metrics,
                    graph_condition=graph_condition,
                    graph_family="eFG" if graph_condition.startswith("eFG") else "STRING",
                    string_threshold="" if graph_condition.startswith("eFG") else graph_condition.split("_", 1)[1],
                    run_group="supplementary_source_comparison",
                    metrics_path=metrics_path,
                    model=args.model,
                    feature_setting=args.feature_setting,
                )
            )

    supplementary_per_run_df = pd.DataFrame(supplementary_rows)
    supplementary_per_run_df.to_csv(summary_dir / "Figure4_supplementary_source_comparison_per_run_metrics.tsv", sep="\t", index=False)
    supplementary_summary_df = aggregate_metrics(
        supplementary_per_run_df,
        ["graph_condition", "graph_family", "model", "feature_setting", "model_feature_setting"],
        METRIC_COLS,
    )
    supplementary_summary_df["error_bar_definition"] = "mean_plus_minus_sd"
    supplementary_summary_df = ordered_source_summary(
        supplementary_summary_df,
        ["STRING_300", "STRING_700", "eFG_HIGH", "eFG_HIGH_MEDIUM", "eFG_ALL"],
    )
    supplementary_summary_df["adapter_summary_path"] = supplementary_summary_df["graph_condition"].map(
        {
            "eFG_HIGH": supplementary_adapter_summaries.get("eFG_HIGH", ""),
            "eFG_HIGH_MEDIUM": supplementary_adapter_summaries.get("eFG_HIGH_MEDIUM", ""),
            "eFG_ALL": supplementary_adapter_summaries.get("eFG_ALL", supplementary_adapter_summaries.get("eFG", "")),
            "STRING_300": "",
            "STRING_700": "",
        }
    )
    supplementary_summary_df.to_csv(summary_dir / "Figure4_supplementary_source_comparison_metrics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
