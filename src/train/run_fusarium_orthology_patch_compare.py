import io
import os
import sys

import numpy as np
import pandas as pd

from src.train.train_support_graph_baseline import (
    require_epgat_env,
    load_yaml,
    set_seed,
    train_one_model,
    build_dgl_graph,
)


SUPPORT_SPECIES = ["human", "scerevisiae", "celegans"]
NEW_FUSARIUM_COLS = [
    "presence_species_count",
    "occupancy",
    "single_copy_fraction",
    "pangenome_core",
    "pangenome_softcore",
    "pangenome_shell",
    "pangenome_cloud",
    "orthology_patch_has_feature",
    "orthology_patch_missing_mask",
    "orthology_patch_ambiguous_mask",
]


def enable_utf8_stdout():
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    except Exception:
        pass


def _read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _load_master():
    return _read_tsv("data_registry/master_label_table.preliminary.tsv")


def _load_wb_alias():
    return _read_tsv("data_registry/wb_alias_mapping.tsv")


def _build_support_patch_features():
    outdir = "outputs/fusarium_feature_patch/support_orthology_patch"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    master = _load_master()
    alias = _load_wb_alias()
    alias["canonical_gene_id"] = alias["WBGene_id"].map(
        dict(
            zip(
                master[master["species"] == "celegans"]["raw_gene_id"],
                master[master["species"] == "celegans"]["canonical_gene_id"],
            )
        )
    ).fillna("")
    alias = alias[alias["canonical_gene_id"].astype(str).ne("")].drop_duplicates(subset=["alias"], keep="first")
    ce_alias_map = dict(zip(alias["alias"], alias["canonical_gene_id"]))
    human_map = dict(
        zip(
            master[master["species"] == "human"]["raw_protein_id"],
            master[master["species"] == "human"]["canonical_gene_id"],
        )
    )
    yeast_map = dict(
        zip(
            master[master["species"] == "scerevisiae"]["raw_protein_id"],
            master[master["species"] == "scerevisiae"]["canonical_gene_id"],
        )
    )

    species_asset = {
        "human": "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/human/Orthologs/orthologs.csv",
        "scerevisiae": "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/yeast/Orthologs/orthologs.csv",
        "celegans": "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/celegans/Orthologs/orthologs.csv",
    }
    species_map = {"human": human_map, "scerevisiae": yeast_map, "celegans": ce_alias_map}

    for species in SUPPORT_SPECIES:
        base = pd.read_csv(
            "outputs/support_graph_features/support_feature_matrix_{}.tsv".format(species),
            sep="\t",
            dtype=str,
        ).fillna("")
        ortho = pd.read_csv(species_asset[species]).fillna(0)
        numeric_cols = [c for c in ortho.columns if c != "Gene"]
        total_species = float(len(numeric_cols))
        ortho["canonical_gene_id"] = ortho["Gene"].astype(str).map(species_map[species]).fillna("")
        ortho = ortho[ortho["canonical_gene_id"].astype(str).ne("")].copy()
        values = ortho[numeric_cols].astype(float)
        presence_species_count = (values > 0).sum(axis=1).astype(float)
        single_copy_count = (values == 1).sum(axis=1).astype(float)
        occupancy = presence_species_count / total_species if total_species > 0 else 0.0
        single_copy_fraction = np.where(presence_species_count > 0, single_copy_count / presence_species_count, 0.0)
        patch = pd.DataFrame(
            {
                "canonical_gene_id": ortho["canonical_gene_id"].astype(str),
                "presence_species_count": presence_species_count,
                "occupancy": occupancy,
                "single_copy_fraction": single_copy_fraction,
            }
        )
        patch["pangenome_core"] = (patch["occupancy"] >= 0.95).astype(float)
        patch["pangenome_softcore"] = ((patch["occupancy"] >= 0.80) & (patch["occupancy"] < 0.95)).astype(float)
        patch["pangenome_shell"] = ((patch["occupancy"] >= 0.15) & (patch["occupancy"] < 0.80)).astype(float)
        patch["pangenome_cloud"] = (patch["occupancy"] < 0.15).astype(float)
        patch["orthology_patch_has_feature"] = 1.0
        patch["orthology_patch_missing_mask"] = 0.0
        patch["orthology_patch_ambiguous_mask"] = 0.0
        patch = patch.groupby("canonical_gene_id", as_index=False).max()
        merged = base.merge(patch, on="canonical_gene_id", how="left")
        for col in NEW_FUSARIUM_COLS:
            merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0.0)
        merged.loc[merged["orthology_patch_has_feature"] <= 0, "orthology_patch_missing_mask"] = 1.0
        merged.to_csv(os.path.join(outdir, "support_feature_matrix_{}.tsv".format(species)), sep="\t", index=False)
    return outdir


def _build_before_after_fusarium_matrices():
    outdir = "outputs/fusarium_feature_patch"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    before = pd.read_csv("outputs/fusarium_graph_inference/fusarium_graph_feature_matrix.tsv", sep="\t", dtype=str).fillna("")
    before.to_csv(os.path.join(outdir, "fusarium_feature_matrix_before_patch.tsv"), sep="\t", index=False)

    orth = _read_tsv("outputs/fusarium_orthology/fusarium_orthology_features.tsv")
    patch = orth[["canonical_gene_id", "presence_species_count", "occupancy", "single_copy_fraction", "pangenome_class", "orthology_join_status"]].copy()
    patch["presence_species_count"] = pd.to_numeric(patch["presence_species_count"], errors="coerce").fillna(0.0)
    patch["occupancy"] = pd.to_numeric(patch["occupancy"], errors="coerce").fillna(0.0)
    patch["single_copy_fraction"] = pd.to_numeric(patch["single_copy_fraction"], errors="coerce").fillna(0.0)
    patch["pangenome_core"] = (patch["pangenome_class"] == "core").astype(float)
    patch["pangenome_softcore"] = (patch["pangenome_class"] == "softcore").astype(float)
    patch["pangenome_shell"] = (patch["pangenome_class"] == "shell").astype(float)
    patch["pangenome_cloud"] = (patch["pangenome_class"] == "cloud").astype(float)
    patch["orthology_patch_has_feature"] = (patch["orthology_join_status"] == "exact").astype(float)
    patch["orthology_patch_missing_mask"] = (patch["orthology_join_status"] == "missing_in_orthofinder_or_unresolved").astype(float)
    patch["orthology_patch_ambiguous_mask"] = (patch["orthology_join_status"] == "ambiguous_multiple_orthogroups").astype(float)
    patch = patch[
        [
            "canonical_gene_id",
            "presence_species_count",
            "occupancy",
            "single_copy_fraction",
            "pangenome_core",
            "pangenome_softcore",
            "pangenome_shell",
            "pangenome_cloud",
            "orthology_patch_has_feature",
            "orthology_patch_missing_mask",
            "orthology_patch_ambiguous_mask",
            "orthology_join_status",
        ]
    ]

    after = before.merge(patch, on="canonical_gene_id", how="left")
    for col in NEW_FUSARIUM_COLS:
        after[col] = pd.to_numeric(after[col], errors="coerce").fillna(0.0)
    after["orthology_join_status"] = after["orthology_join_status"].replace("", pd.NA).fillna("missing_in_orthofinder_or_unresolved")
    after.to_csv(os.path.join(outdir, "fusarium_feature_matrix_after_patch.tsv"), sep="\t", index=False)

    val = {
        "row_count_before": len(before),
        "row_count_after": len(after),
        "join_success_count": int((after["orthology_patch_has_feature"] > 0).sum()),
        "unresolved_count": int((after["orthology_patch_missing_mask"] > 0).sum()),
        "ambiguous_count": int((after["orthology_patch_ambiguous_mask"] > 0).sum()),
        "row_count_match": len(before) == len(after),
    }
    validation_rows = []
    for col in ["presence_species_count", "occupancy", "single_copy_fraction", "pangenome_core", "pangenome_softcore", "pangenome_shell", "pangenome_cloud"]:
        missing_rate = float((pd.to_numeric(after[col], errors="coerce").fillna(0.0) == 0).mean())
        validation_rows.append(
            {
                "feature_name": col,
                "missing_rate": missing_rate,
                "join_success_count": val["join_success_count"],
                "unresolved_count": val["unresolved_count"],
                "ambiguous_count": val["ambiguous_count"],
                "encoding": "numeric" if not col.startswith("pangenome_") else "one_hot",
            }
        )
    pd.DataFrame(validation_rows).to_csv(os.path.join(outdir, "fusarium_feature_patch_validation_audit.tsv"), sep="\t", index=False)
    schema_rows = [
        {"feature_name": "presence_species_count", "encoding": "numeric", "source": "fusarium_orthology_features.tsv"},
        {"feature_name": "occupancy", "encoding": "numeric", "source": "fusarium_orthology_features.tsv"},
        {"feature_name": "single_copy_fraction", "encoding": "numeric", "source": "fusarium_orthology_features.tsv"},
        {"feature_name": "pangenome_core", "encoding": "one_hot", "source": "pangenome_class core"},
        {"feature_name": "pangenome_softcore", "encoding": "one_hot", "source": "pangenome_class softcore"},
        {"feature_name": "pangenome_shell", "encoding": "one_hot", "source": "pangenome_class shell"},
        {"feature_name": "pangenome_cloud", "encoding": "one_hot", "source": "pangenome_class cloud"},
        {"feature_name": "orthology_patch_has_feature", "encoding": "mask", "source": "orthology_join_status"},
        {"feature_name": "orthology_patch_missing_mask", "encoding": "mask", "source": "orthology_join_status"},
        {"feature_name": "orthology_patch_ambiguous_mask", "encoding": "mask", "source": "orthology_join_status"},
    ]
    pd.DataFrame(schema_rows).to_csv(os.path.join(outdir, "fusarium_feature_patch_schema.tsv"), sep="\t", index=False)
    md_lines = [
        "# Fusarium Feature Patch Validation",
        "",
        "- row_count_before = {}".format(val["row_count_before"]),
        "- row_count_after = {}".format(val["row_count_after"]),
        "- join_success_count = {}".format(val["join_success_count"]),
        "- unresolved_count = {}".format(val["unresolved_count"]),
        "- ambiguous_count = {}".format(val["ambiguous_count"]),
        "- pangenome_class encoding = one-hot(core, softcore, shell, cloud)",
    ]
    open(os.path.join(outdir, "fusarium_feature_patch_validation_audit.md"), "w", encoding="utf-8").write("\n".join(md_lines))
    return outdir


def _load_fusarium_graph_features(feature_path, feature_columns):
    df = pd.read_csv(feature_path, sep="\t", dtype=str).fillna("")
    numeric_cols = []
    for col in feature_columns:
        if col == "embedding_vector":
            numeric_cols.extend([c for c in df.columns if c.startswith("emb_")])
        elif col in df.columns:
            numeric_cols.append(col)
    x = df[numeric_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)
    return df, x, numeric_cols


def _infer_scores(model, feature_path, feature_columns):
    import torch
    edges = pd.read_csv("outputs/fusarium_graph_inference/fusarium_edges_for_inference.tsv", sep="\t", dtype=str).fillna("")
    nodes_df, x_np, cols = _load_fusarium_graph_features(feature_path, feature_columns)
    node_index = dict(zip(nodes_df["canonical_gene_id"].astype(str), range(len(nodes_df))))
    edges = edges.copy()
    edges["source_node_index"] = edges["source_canonical_gene_id"].map(node_index)
    edges["target_node_index"] = edges["target_canonical_gene_id"].map(node_index)
    edges = edges.dropna(subset=["source_node_index", "target_node_index"]).copy()
    graph = build_dgl_graph(edges, len(nodes_df))
    feats = torch.tensor(x_np, dtype=torch.float32)
    model.eval()
    with torch.no_grad():
        logits = model(graph, feats).squeeze()
        scores = torch.sigmoid(logits).cpu().numpy()
    out = nodes_df[["canonical_gene_id"]].copy()
    out["prediction_score"] = scores
    out = out.sort_values("prediction_score", ascending=False, kind="stable").reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def _ranking_metrics(score_df):
    benchmark_files = {
        "broad79": "outputs/baseline_dataset/fgraminearum_broad79.tsv",
        "strict29": "outputs/baseline_dataset/fgraminearum_strict29.tsv",
        "conflict8": "outputs/baseline_dataset/fgraminearum_conflict8.tsv",
    }
    rows = []
    for subset, path in benchmark_files.items():
        positives = set(_read_tsv(path)["canonical_gene_id"].astype(str))
        ranked = score_df["canonical_gene_id"].astype(str).tolist()
        for k in [50, 100, 200, 500]:
            topk = ranked[:k]
            hits = [g for g in topk if g in positives]
            hit_count = len(hits)
            rows.append(
                {
                    "subset": subset,
                    "k": k,
                    "hit_count": hit_count,
                    "precision_at_k": hit_count / float(k),
                    "recall_at_k": hit_count / float(len(positives)) if positives else 0.0,
                }
            )
    return pd.DataFrame(rows)


def main():
    enable_utf8_stdout()
    require_epgat_env()
    print("【开始 Fusarium orthology feature patch 接入与排名比较】")
    print("当前固定模型：")
    print("  • GraphSAGE")
    print("  • hidden_dim = 64")
    print("  • num_layers = 2")
    print("  • dropout = 0.2")
    print("新增特征：")
    print("  • presence_species_count")
    print("  • occupancy")
    print("  • single_copy_fraction")
    print("  • pangenome_class")
    print("当前规则：")
    print("  • pangenome_class 完全复用旧脚本定义")
    print("  • 不重定义阈值")

    patch_outdir = _build_before_after_fusarium_matrices()
    support_patch_dir = _build_support_patch_features()

    print("【开始构建 patch 前 feature matrix】")
    print("【开始构建 patch 后 feature matrix】")

    base_cfg = load_yaml("configs/support_graph_baseline.yaml")
    fixed_cfg = dict(base_cfg)
    fixed_cfg["species_sets"] = {"default": SUPPORT_SPECIES}
    fixed_cfg["species_loss_weights"] = {"human": 1.0, "scerevisiae": 1.0, "celegans": 1.0}
    fixed_cfg["hidden_dim"] = 64
    fixed_cfg["num_layers"] = 2
    fixed_cfg["dropout"] = 0.2
    fixed_cfg["feature_normalization"] = False
    fixed_cfg["prior_feature_dir"] = "outputs/support_prior"

    print("【开始运行 patch 前 Fusarium 排名】")
    fixed_cfg["orthology_feature_dir"] = "outputs/support_graph_features"
    fixed_cfg["feature_scope"] = {
        "embedding": True,
        "expression": False,
        "orthology": True,
        "localization": False,
        "prior_score": True,
    }
    set_seed(int(fixed_cfg.get("seed", 20260403)))
    before_train = train_one_model("GraphSAGE", fixed_cfg, "default")
    before_scores = _infer_scores(
        before_train["model_object"],
        os.path.join(patch_outdir, "fusarium_feature_matrix_before_patch.tsv"),
        before_train["feature_columns"],
    )
    before_scores.to_csv(os.path.join(patch_outdir, "graphsage_before_patch_scores.tsv"), sep="\t", index=False)
    before_metrics = _ranking_metrics(before_scores)
    before_metrics.to_csv(os.path.join(patch_outdir, "graphsage_before_patch_metrics.tsv"), sep="\t", index=False)

    print("【开始运行 patch 后 Fusarium 排名】")
    fixed_cfg["orthology_feature_dir"] = support_patch_dir
    set_seed(int(fixed_cfg.get("seed", 20260403)))
    after_train = train_one_model("GraphSAGE", fixed_cfg, "default")
    after_scores = _infer_scores(
        after_train["model_object"],
        os.path.join(patch_outdir, "fusarium_feature_matrix_after_patch.tsv"),
        after_train["feature_columns"],
    )
    after_scores.to_csv(os.path.join(patch_outdir, "graphsage_after_patch_scores.tsv"), sep="\t", index=False)
    after_metrics = _ranking_metrics(after_scores)
    after_metrics.to_csv(os.path.join(patch_outdir, "graphsage_after_patch_metrics.tsv"), sep="\t", index=False)

    print("【开始计算 strict29 / broad79 / conflict8 top-k 命中】")
    comp = before_metrics.merge(after_metrics, on=["subset", "k"], suffixes=("_before", "_after"))
    comp["delta_hit_count"] = comp["hit_count_after"] - comp["hit_count_before"]
    comp["delta_precision_at_k"] = comp["precision_at_k_after"] - comp["precision_at_k_before"]
    comp["delta_recall_at_k"] = comp["recall_at_k_after"] - comp["recall_at_k_before"]
    comp.to_csv(os.path.join(patch_outdir, "graphsage_patch_comparison.tsv"), sep="\t", index=False)

    summary_lines = [
        "# GraphSAGE Patch Comparison",
        "",
        "- strict29 top100 before = {}".format(int(comp[(comp["subset"] == "strict29") & (comp["k"] == 100)]["hit_count_before"].iloc[0])),
        "- strict29 top100 after = {}".format(int(comp[(comp["subset"] == "strict29") & (comp["k"] == 100)]["hit_count_after"].iloc[0])),
        "- broad79 top100 before = {}".format(int(comp[(comp["subset"] == "broad79") & (comp["k"] == 100)]["hit_count_before"].iloc[0])),
        "- broad79 top100 after = {}".format(int(comp[(comp["subset"] == "broad79") & (comp["k"] == 100)]["hit_count_after"].iloc[0])),
        "- conflict8 top100 before = {}".format(int(comp[(comp["subset"] == "conflict8") & (comp["k"] == 100)]["hit_count_before"].iloc[0])),
        "- conflict8 top100 after = {}".format(int(comp[(comp["subset"] == "conflict8") & (comp["k"] == 100)]["hit_count_after"].iloc[0])),
    ]
    open(os.path.join(patch_outdir, "graphsage_patch_comparison.md"), "w", encoding="utf-8").write("\n".join(summary_lines))
    open("101_fusarium_orthology_patch_results.md", "w", encoding="utf-8").write("\n".join(summary_lines))
    next_lines = [
        "# 102 Next Step After Fusarium Orthology Patch",
        "",
        "- 如果 patch 改善了 strict29 / broad79 顶部命中，则应将其固定为默认 Fusarium-side GraphSAGE 特征集。",
        "- 然后进入真实 Fusarium candidate ranking 主结果阶段。",
    ]
    open("102_next_step_after_fusarium_orthology_patch.md", "w", encoding="utf-8").write("\n".join(next_lines))

    print("【Fusarium orthology patch 比较完成】")
    print("结果摘要：")
    for subset in ["strict29", "broad79"]:
        row = comp[(comp["subset"] == subset) & (comp["k"] == 100)].iloc[0]
        print("  • patch 前 {} top100 命中 = {}".format(subset, int(row["hit_count_before"])))
        print("  • patch 后 {} top100 命中 = {}".format(subset, int(row["hit_count_after"])))
    improved = (comp["delta_recall_at_k"] > 0).any() or (comp["delta_precision_at_k"] > 0).any()
    print("当前结论：")
    print("  • 是否建议将 Fusarium orthology patch 纳入默认特征集：{}".format("是" if improved else "否"))


if __name__ == "__main__":
    main()
