from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml
import numpy as np

from src.registry.support_graph_registry import SupportGraphRegistryRecord


def load_config(path: str | Path) -> dict:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def main() -> None:
    config = load_config("configs/support_graph_build.yaml")
    out = Path("outputs/support_graphs")
    out.mkdir(parents=True, exist_ok=True)

    admission = _read_tsv("outputs/support_graphs/support_graph_admission.tsv")
    support_samples = _read_tsv("outputs/baseline_dataset/support_supervised_samples.tsv")
    embedding_manifest = _read_tsv("outputs/baseline_dataset/embedding_manifest.pooled.tsv")
    valid_embeddings = embedding_manifest[
        embedding_manifest["exists"].astype(str).str.lower().isin({"true", "1", "yes"})
        & ~embedding_manifest["needs_manual_review"].astype(str).str.lower().isin({"true", "1", "yes"})
    ].copy()

    registry_rows = []
    summary_lines = [
        "# Support Graph Build Summary",
        "",
        "## Direct Answers",
    ]

    for _, row in admission.iterrows():
        species = row["species"]
        print(f"[support_graph_build] species={species} start")
        use_in_training = row["is_support_graph_usable"].lower() in {"true", "1", "yes"}
        graph_class = row["graph_completeness_class"]
        default_edge_weight = (
            float(config["graph_weight_rules"]["full"])
            if graph_class == "graph-complete-support"
            else float(config["graph_weight_rules"]["partial"])
        )
        edge_usage_mode = "full" if graph_class == "graph-complete-support" else "downweighted"

        species_samples = support_samples[support_samples["species"] == species].copy()
        species_embeddings = valid_embeddings[valid_embeddings["species"] == species].copy()
        embedding_ids = set(species_embeddings["canonical_gene_id"].astype(str))

        edge_input = _read_tsv(f"outputs/support_graphs/{species}_ppi_edges_bridged.tsv")
        print(f"[support_graph_build] species={species} raw_edge_rows={len(edge_input)}")
        edge_input["source_canonical_gene_id"] = edge_input["source_canonical_gene_id"].astype(str)
        edge_input["target_canonical_gene_id"] = edge_input["target_canonical_gene_id"].astype(str)
        edge_input = edge_input[
            edge_input["source_canonical_gene_id"].ne("")
            & edge_input["target_canonical_gene_id"].ne("")
        ].copy()
        print(f"[support_graph_build] species={species} fully_bridged_edge_rows={len(edge_input)}")
        edge_input = edge_input[
            edge_input["source_canonical_gene_id"].isin(embedding_ids)
            & edge_input["target_canonical_gene_id"].isin(embedding_ids)
        ].copy()
        edge_input = edge_input[
            edge_input["source_canonical_gene_id"].ne(edge_input["target_canonical_gene_id"])
        ].copy()
        if not edge_input.empty:
            edge_input["source_gene_id"] = edge_input["source_canonical_gene_id"]
            edge_input["target_gene_id"] = edge_input["target_canonical_gene_id"]
            src_vals = edge_input["source_gene_id"].to_numpy(dtype=object)
            dst_vals = edge_input["target_gene_id"].to_numpy(dtype=object)
            ordered_src = np.where(src_vals <= dst_vals, src_vals, dst_vals)
            ordered_dst = np.where(src_vals <= dst_vals, dst_vals, src_vals)
            edge_input["source_gene_id"] = ordered_src
            edge_input["target_gene_id"] = ordered_dst
            edge_input["raw_edge_weight"] = edge_input["edge_weight"].astype(float)
            edge_input = edge_input.sort_values(
                ["source_gene_id", "target_gene_id", "raw_edge_weight"],
                ascending=[True, True, False],
                kind="stable",
            )
            edge_input = edge_input.drop_duplicates(
                subset=["source_gene_id", "target_gene_id"], keep="first"
            )[
                ["source_gene_id", "target_gene_id", "edge_weight", "raw_edge_weight"]
            ].reset_index(drop=True)
            print(f"[support_graph_build] species={species} dedup_edge_rows={len(edge_input)}")
        else:
            edge_input = pd.DataFrame(columns=["source_gene_id", "target_gene_id", "edge_weight", "raw_edge_weight"])

        graph_nodes = sorted(set(edge_input["source_gene_id"].astype(str)) | set(edge_input["target_gene_id"].astype(str)))
        if not graph_nodes:
            graph_nodes = sorted(set(species_embeddings["canonical_gene_id"].astype(str)))

        nodes_df = pd.DataFrame({"canonical_gene_id": graph_nodes})
        nodes_df["species"] = species
        nodes_df["graph_node_id"] = nodes_df["canonical_gene_id"]
        nodes_df["is_in_embedding_universe"] = True
        nodes_df["graph_completeness_class"] = graph_class
        nodes_df["node_mapping_rate"] = row["safe_bridge_ratio"]
        nodes_df["edge_retention_rate"] = row["edge_retention_rate"]
        nodes_df["graph_weight_class"] = "full" if graph_class == "graph-complete-support" else "partial"
        nodes_df["use_in_training"] = use_in_training
        nodes_df = nodes_df[
            [
                "species",
                "canonical_gene_id",
                "graph_node_id",
                "is_in_embedding_universe",
                "graph_completeness_class",
                "node_mapping_rate",
                "edge_retention_rate",
                "graph_weight_class",
                "use_in_training",
            ]
        ]

        edge_df = edge_input.copy()
        edge_df["species"] = species
        edge_df["graph_completeness_class"] = graph_class
        edge_df["edge_usage_mode"] = edge_usage_mode
        if "raw_edge_weight" not in edge_df.columns:
            edge_df["raw_edge_weight"] = edge_df["edge_weight"].astype(float) if not edge_df.empty else []
        edge_df["edge_weight"] = default_edge_weight
        edge_df = edge_df[
            [
                "species",
                "source_gene_id",
                "target_gene_id",
                "edge_weight",
                "graph_completeness_class",
                "edge_usage_mode",
                "raw_edge_weight",
            ]
        ]

        feature_df = species_embeddings[species_embeddings["canonical_gene_id"].isin(graph_nodes)].copy()
        feature_df = feature_df[["species", "canonical_gene_id", "embedding_source"]].copy()
        feature_df["has_embedding"] = True
        feature_df["has_expression_feature"] = False
        feature_df["has_orthology_feature"] = False
        feature_df["has_localization_feature"] = False
        feature_df["has_prior_score"] = False
        feature_df["feature_ready"] = True
        feature_df = feature_df[
            [
                "species",
                "canonical_gene_id",
                "has_embedding",
                "embedding_source",
                "has_expression_feature",
                "has_orthology_feature",
                "has_localization_feature",
                "has_prior_score",
                "feature_ready",
            ]
        ]

        label_df = species_samples[species_samples["canonical_gene_id"].isin(graph_nodes)].copy()
        label_df["is_labeled"] = label_df["gold_label"].astype(str).isin({"0", "1"}) & label_df["label_status"].astype(str).eq("gold")
        label_df["label_source"] = label_df["notes"]
        label_df["label_value"] = label_df["gold_label"]
        label_df["usable_for_support_supervision"] = label_df["is_labeled"]
        label_df = label_df[
            [
                "species",
                "canonical_gene_id",
                "is_labeled",
                "label_source",
                "label_value",
                "usable_for_support_supervision",
            ]
        ]

        node_path = out / f"{species}_nodes.tsv"
        edge_path = out / f"{species}_edges_for_training.tsv"
        feature_path = out / f"{species}_node_feature_manifest.tsv"
        label_path = out / f"{species}_label_manifest.tsv"

        nodes_df.to_csv(node_path, sep="\t", index=False)
        edge_df.to_csv(edge_path, sep="\t", index=False)
        feature_df.to_csv(feature_path, sep="\t", index=False)
        label_df.to_csv(label_path, sep="\t", index=False)
        print(
            f"[support_graph_build] species={species} wrote nodes={len(nodes_df)} edges={len(edge_df)} "
            f"features={len(feature_df)} labels={len(label_df)}"
        )

        note = (
            "partial-support species: include with downweighted edges"
            if graph_class == "partial-support-graph"
            else "graph-complete support species"
        )
        registry_rows.append(
            SupportGraphRegistryRecord(
                species=species,
                graph_completeness_class=graph_class,
                node_table_path=str(node_path),
                edge_table_path=str(edge_path),
                node_feature_manifest_path=str(feature_path),
                label_manifest_path=str(label_path),
                edge_weight_default=default_edge_weight,
                use_in_graph_training=bool(use_in_training),
                notes=note,
            ).to_dict()
        )

        summary_lines.append(
            f"- {species}: nodes={len(nodes_df)}, edges={len(edge_df)}, class={graph_class}, "
            f"edge_weight_default={default_edge_weight}, use_in_training={str(use_in_training).lower()}"
        )

    registry_df = pd.DataFrame(registry_rows)
    registry_df.to_csv(out / "support_graph_registry.tsv", sep="\t", index=False)

    summary_lines.extend(
        [
            "",
            "## Interpretation",
            "- human 与 scerevisiae 作为 graph-complete-support 纳入训练输入。",
            "- celegans 作为 partial-support-graph 纳入训练输入，并保留降权 caveat。",
            "- 当前图输入构建已完成；真正缺失的是下游 GNN / GAT / GraphSAGE 训练实现，而不是 graph-ready 结构化输入。",
        ]
    )
    (Path("68_support_graph_build_summary.md")).write_text("\n".join(summary_lines), encoding="utf-8")

    next_lines = [
        "# 69 Next Step After Support Graph Build",
        "",
        "下一步应进入 graph-aware support-species model input assembly 与首个 GNN / GraphSAGE baseline training。",
        "",
        "建议起点：",
        "- 先用 embeddings + PPI 构建最小 graph baseline",
        "- human 与 scerevisiae 作为 full graph support",
        "- celegans 作为 partial support graph 保留降权输入",
        "- expression / orthology / localization 后续再逐步接入，不在本阶段伪造为 ready",
    ]
    Path("69_next_step_after_support_graph_build.md").write_text("\n".join(next_lines), encoding="utf-8")

    Path("outputs/support_graphs/support_graph_build_config_snapshot.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False),
        encoding="utf-8",
    )

    print("【support graph 构建完成】")
    print("已纳入物种：")
    for _, row in registry_df.iterrows():
        label = "完整图支持（graph-complete-support）" if row["graph_completeness_class"] == "graph-complete-support" else "部分图支持（partial-support-graph，训练时降权）"
        print(f"  - {row['species']}：{label}")
    print("已生成输出：")
    print("  - 每个物种的节点表")
    print("  - 每个物种的训练边表")
    print("  - 每个物种的特征清单")
    print("  - 每个物种的标签清单")
    print("  - support graph 总注册表")
    print("当前结论：")
    print("  - 已完成 graph-ready 输入构建")
    print("  - 下一步可进入 GNN / GAT / GraphSAGE 训练阶段")


if __name__ == "__main__":
    main()
