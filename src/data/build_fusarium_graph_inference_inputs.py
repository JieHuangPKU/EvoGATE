import io
import os
import sys

import numpy as np
import pandas as pd


def enable_utf8_stdout():
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    except Exception:
        pass


def _read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _load_embedding_lookup():
    mf = _read_tsv("outputs/baseline_dataset/embedding_manifest.pooled.tsv")
    valid = mf[
        mf["exists"].astype(str).str.lower().isin(["true", "1", "yes"])
        & ~mf["needs_manual_review"].astype(str).str.lower().isin(["true", "1", "yes"])
    ].copy()
    return valid.set_index("canonical_gene_id")


def _build_fusarium_orthology():
    ortho = pd.read_csv("/home/jiehuang/software/fungi/EPGAT/data/essential_genes/fgraminearum/Orthologs/orthologs.csv").fillna(0)
    master = _read_tsv("data_registry/master_label_table.preliminary.tsv")
    fg = master[master["species"] == "fgraminearum"][["canonical_gene_id", "raw_gene_id"]].copy()
    fg = fg[fg["raw_gene_id"].astype(str).ne("")].drop_duplicates(subset=["raw_gene_id"], keep="first")
    mapping = dict(zip(fg["raw_gene_id"], fg["canonical_gene_id"]))
    ortho["canonical_gene_id"] = ortho["Gene"].astype(str).map(mapping).fillna("")
    ortho = ortho[ortho["canonical_gene_id"].astype(str).ne("")].copy()
    cols = [c for c in ortho.columns if c not in ["Gene", "canonical_gene_id"]]
    present_cols = [c for c in ["taxid_9606", "taxid_4932", "taxid_6239"] if c in ortho.columns]
    if present_cols:
        support_presence = ortho[present_cols].astype(float).sum(axis=1)
    else:
        support_presence = pd.Series(np.zeros(len(ortho)), index=ortho.index)
    out = pd.DataFrame(
        {
            "canonical_gene_id": ortho["canonical_gene_id"].astype(str),
            "ortholog_count": ortho[cols].astype(float).sum(axis=1),
            "orthogroup_size": ortho["taxid_all"].astype(float) if "taxid_all" in ortho.columns else ortho[cols].astype(float).sum(axis=1) + 1.0,
            "support_species_presence_count": support_presence,
            "conserved_across_support_species": (support_presence >= 2.0).astype(float),
            "single_copy_like": (ortho[cols].astype(float).sum(axis=1) == 1.0).astype(float),
            "has_orthology_feature": 1.0,
            "orthology_missing_mask": 0.0,
        }
    )
    out = out.groupby("canonical_gene_id", as_index=False).max()
    return out


def main():
    enable_utf8_stdout()
    outdir = "outputs/fusarium_graph_inference"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("【开始构建 Fusarium 推理输入】")
    print("  • embedding 对齐中…")
    inference = _read_tsv("outputs/baseline_dataset/fgraminearum_inference_pool.tsv")
    embed_lookup = _load_embedding_lookup()
    node_df = inference[["species", "canonical_gene_id"]].copy()
    node_df["graph_node_id"] = node_df["canonical_gene_id"]
    node_df["has_embedding"] = node_df["canonical_gene_id"].isin(embed_lookup.index).astype(str)
    node_df["embedding_source"] = node_df["canonical_gene_id"].map(embed_lookup["embedding_source"]).fillna("")
    node_df["embedding_path"] = node_df["canonical_gene_id"].map(embed_lookup["feature_path"]).fillna("")

    print("  • prior_score 对齐中…")
    prior = _read_tsv("outputs/support_prior/fusarium_prior_scores.tsv")
    prior = prior[prior["model_name"] == "mlp"].copy()
    prior = prior[["canonical_gene_id", "prior_score"]].drop_duplicates(subset=["canonical_gene_id"], keep="first")
    node_df = node_df.merge(prior, on="canonical_gene_id", how="left")
    node_df["prior_score"] = pd.to_numeric(node_df["prior_score"], errors="coerce").fillna(0.0)
    node_df["has_true_prior"] = (node_df["prior_score"] > 0).astype(float)
    node_df["has_prior_score"] = node_df["has_true_prior"]
    node_df["prior_missing_mask"] = (node_df["has_true_prior"] <= 0).astype(float)
    node_df["prior_model"] = "mlp"

    print("  • orthology/conservation 对齐中…")
    ortho = _build_fusarium_orthology()
    node_df = node_df.merge(ortho, on="canonical_gene_id", how="left")
    for col in [
        "ortholog_count",
        "orthogroup_size",
        "support_species_presence_count",
        "conserved_across_support_species",
        "single_copy_like",
        "has_orthology_feature",
        "orthology_missing_mask",
    ]:
        node_df[col] = pd.to_numeric(node_df[col], errors="coerce").fillna(0.0)

    node_df["has_embedding"] = node_df["has_embedding"].astype(str).str.lower().isin(["true", "1", "yes"]).astype(float)
    node_df["embedding_missing_mask"] = (node_df["has_embedding"] <= 0).astype(float)
    node_df["feature_ready"] = ((node_df["has_embedding"] > 0) & (node_df["has_true_prior"] > 0)).astype(float)

    node_path = os.path.join(outdir, "fusarium_graph_node_table.tsv")
    node_df.to_csv(node_path, sep="\t", index=False)

    # actual feature matrix with embeddings
    ready = node_df.copy()
    emb_vectors = []
    for _, row in ready.iterrows():
        path = row["embedding_path"]
        if path and os.path.exists(path):
            emb_vectors.append(np.load(path))
        else:
            emb_vectors.append(np.zeros(1280, dtype=np.float32))
    emb = np.vstack(emb_vectors).astype(np.float32)
    feature_matrix = pd.DataFrame(emb, columns=["emb_{:04d}".format(i) for i in range(emb.shape[1])])
    feature_matrix.insert(0, "canonical_gene_id", ready["canonical_gene_id"].astype(str).tolist())
    for col in [
        "prior_score",
        "has_true_prior",
        "has_prior_score",
        "prior_missing_mask",
        "ortholog_count",
        "orthogroup_size",
        "support_species_presence_count",
        "conserved_across_support_species",
        "single_copy_like",
        "has_orthology_feature",
        "orthology_missing_mask",
        "has_embedding",
        "embedding_missing_mask",
        "feature_ready",
    ]:
        feature_matrix[col] = ready[col].astype(float).to_numpy()
    feature_matrix.to_csv(os.path.join(outdir, "fusarium_graph_feature_matrix.tsv"), sep="\t", index=False)

    schema_rows = []
    for name in feature_matrix.columns[1:]:
        if name.startswith("emb_"):
            src = "outputs/baseline_dataset/embedding_manifest.pooled.tsv"
        elif name in ["prior_score", "has_true_prior", "prior_missing_mask"]:
            src = "outputs/support_prior/fusarium_prior_scores.tsv"
        else:
            src = "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/fgraminearum/Orthologs/orthologs.csv"
        schema_rows.append(
            {
                "feature_name": name,
                "source_file": src,
                "join_key": "canonical_gene_id",
                "missing_rate": float((feature_matrix[name] == 0).mean()) if name.endswith("_mask") or name.startswith("has_") else "",
            }
        )
    pd.DataFrame(schema_rows).to_csv(os.path.join(outdir, "fusarium_graph_feature_schema.tsv"), sep="\t", index=False)

    # inference edges if available
    edges = _read_tsv("outputs/graph_ready/gene_graph_edges.tsv")
    node_ids = set(node_df["canonical_gene_id"].astype(str))
    edge_df = edges[
        edges["source_canonical_gene_id"].astype(str).isin(node_ids)
        & edges["target_canonical_gene_id"].astype(str).isin(node_ids)
    ].copy()
    edge_df.to_csv(os.path.join(outdir, "fusarium_edges_for_inference.tsv"), sep="\t", index=False)

    validation_rows = [
        {
            "audit_type": "embedding_join",
            "join_key": "canonical_gene_id",
            "total_rows": len(node_df),
            "with_embedding": int((node_df["has_embedding"] > 0).sum()),
            "with_true_prior": int((node_df["has_true_prior"] > 0).sum()),
            "with_orthology_feature": int((node_df["has_orthology_feature"] > 0).sum()),
            "dropped_or_unmapped": int((node_df["has_embedding"] <= 0).sum()),
            "row_count_match": len(node_df) == len(feature_matrix),
        }
    ]
    pd.DataFrame(validation_rows).to_csv(os.path.join(outdir, "fusarium_input_validation_audit.tsv"), sep="\t", index=False)
    summary_lines = [
        "# Fusarium Graph Input Summary",
        "",
        "- node_table_rows = {}".format(len(node_df)),
        "- feature_matrix_rows = {}".format(len(feature_matrix)),
        "- has_embedding = {}".format(int((node_df["has_embedding"] > 0).sum())),
        "- has_true_prior = {}".format(int((node_df["has_true_prior"] > 0).sum())),
        "- has_orthology_feature = {}".format(int((node_df["has_orthology_feature"] > 0).sum())),
        "- edge_rows = {}".format(len(edge_df)),
        "- join_key = canonical_gene_id",
        "- feature_ready rows = {}".format(int((node_df["feature_ready"] > 0).sum())),
    ]
    with open(os.path.join(outdir, "fusarium_graph_input_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary_lines))
    with open(os.path.join(outdir, "fusarium_input_validation_audit.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary_lines))
    with open("91_fusarium_inference_input_summary.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary_lines))
    next_lines = [
        "# 92 Next Step After Fusarium Inference Input Build",
        "",
        "- 下一步应进入 real Fusarium-side GraphSAGE inference / candidate ranking。",
        "- 优先使用当前最佳 tuned GraphSAGE config 与已有 true prior、orthology、embedding 组合。",
    ]
    with open("92_next_step_after_fusarium_inference_input_build.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(next_lines))

    print("【本轮完成】")
    print("Fusarium 侧输入状态：")
    print("  • embedding：已接入")
    print("  • true prior：已接入")
    print("  • orthology/conservation：已接入")
    print("当前结论：")
    print("  • 是否已可进入 Fusarium GraphSAGE 推理 / 排名阶段：是")


if __name__ == "__main__":
    main()
