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


def _safe_prior_from_edges(species):
    path = "outputs/support_graphs/{}_edges_for_training.tsv".format(species)
    edges = _read_tsv(path)
    nodes = _read_tsv("outputs/support_graphs/{}_nodes.tsv".format(species))
    if "raw_edge_weight" not in edges.columns:
        raise ValueError("边表缺少 raw_edge_weight: {}".format(path))

    work = edges[["source_gene_id", "target_gene_id", "raw_edge_weight"]].copy()
    work["raw_edge_weight"] = pd.to_numeric(work["raw_edge_weight"], errors="coerce").fillna(0.0)
    src = work.rename(columns={"source_gene_id": "canonical_gene_id"})[["canonical_gene_id", "raw_edge_weight"]]
    dst = work.rename(columns={"target_gene_id": "canonical_gene_id"})[["canonical_gene_id", "raw_edge_weight"]]
    both = pd.concat([src, dst], ignore_index=True)
    agg = both.groupby("canonical_gene_id", as_index=False)["raw_edge_weight"].sum()
    agg = agg.rename(columns={"raw_edge_weight": "weighted_degree_sum"})
    merged = nodes[["species", "canonical_gene_id"]].merge(agg, on="canonical_gene_id", how="left")
    merged["weighted_degree_sum"] = pd.to_numeric(merged["weighted_degree_sum"], errors="coerce").fillna(0.0)
    merged["has_prior_score"] = (merged["weighted_degree_sum"] > 0).astype(float)
    transformed = np.log1p(merged["weighted_degree_sum"].to_numpy(dtype=np.float64))
    if len(transformed) and float(transformed.max()) > float(transformed.min()):
        norm = (transformed - transformed.min()) / (transformed.max() - transformed.min())
    else:
        norm = np.zeros_like(transformed, dtype=np.float64)
    merged["prior_score"] = norm.astype(np.float32)
    merged["prior_missing_mask"] = (merged["has_prior_score"] <= 0).astype(float)
    return merged[["species", "canonical_gene_id", "prior_score", "has_prior_score", "prior_missing_mask"]]


def main():
    enable_utf8_stdout()
    outdir = "outputs/support_graph_prior"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    rows = []
    species_map = {
        "human": "human",
        "scerevisiae": "yeast",
        "celegans": "celegans",
    }

    for species, epgat_species in species_map.items():
        ogee = "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/{}/EssentialGenes/ogee.csv".format(epgat_species)
        rows.append(
            {
                "species": species,
                "file_path": ogee,
                "source_type": "legacy_label_like_asset",
                "candidate_score_column": "Label",
                "join_key": "Gene",
                "directly_usable": "false",
                "needs_remapping": "yes",
                "safe_as_prior": "false",
                "notes": "OGEE is label-like and would leak supervision if used as prior_score.",
            }
        )
        edge_path = "outputs/support_graphs/{}_edges_for_training.tsv".format(species)
        rows.append(
            {
                "species": species,
                "file_path": edge_path,
                "source_type": "derived_proxy_prior_from_support_graph",
                "candidate_score_column": "raw_edge_weight",
                "join_key": "source_gene_id / target_gene_id -> canonical_gene_id",
                "directly_usable": "false",
                "needs_remapping": "no",
                "safe_as_prior": "true",
                "notes": "Safe numeric proxy prior derived without labels: log1p(weighted degree sum) then min-max normalized per species.",
            }
        )

    audit = pd.DataFrame(rows)
    audit.to_csv(os.path.join(outdir, "prior_asset_audit.tsv"), sep="\t", index=False)

    lines = [
        "# Support Prior Asset Audit",
        "",
        "## 结论",
        "- 当前没有可直接安全复用的 numeric prior_score 文件。",
        "- OGEE/EssentialGenes 资产是 label-like，不可直接作为 prior_score。",
        "- 本轮采用的安全 prior 是从当前 support graph 边表无监督导出的 weighted-degree proxy prior。",
        "- 具体变换：per-node sum(raw_edge_weight) -> log1p -> per-species min-max normalization。",
    ]
    with open(os.path.join(outdir, "prior_asset_audit.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    with open("83_support_prior_asset_audit.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    for species in ["human", "scerevisiae", "celegans"]:
        prior = _safe_prior_from_edges(species)
        prior.to_csv(os.path.join(outdir, "support_prior_matrix_{}.tsv".format(species)), sep="\t", index=False)

    print("【开始 support prior 特征审计】")
    print("当前默认物种：")
    print("  • human")
    print("  • scerevisiae")
    print("  • celegans")
    print("当前固定模型：")
    print("  • GraphSAGE")
    print("当前固定特征：")
    print("  • embedding")
    print("  • orthology/conservation")
    print("  • missing masks")
    print("本轮目标：")
    print("  • 接入 prior_score 并比较增益")


if __name__ == "__main__":
    main()
