import io
import os
import sys

import pandas as pd


def enable_utf8_stdout():
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    except Exception:
        pass


def _read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _load_master_label_table():
    path = "data_registry/master_label_table.preliminary.tsv"
    return _read_tsv(path)


def _load_wb_alias():
    path = "data_registry/wb_alias_mapping.tsv"
    df = _read_tsv(path)
    return df


def _human_raw_protein_to_canonical(master_df):
    sdf = master_df[master_df["species"] == "human"].copy()
    sdf = sdf[sdf["raw_protein_id"].astype(str).ne("")]
    return dict(zip(sdf["raw_protein_id"], sdf["canonical_gene_id"]))


def _yeast_raw_protein_to_canonical(master_df):
    sdf = master_df[master_df["species"] == "scerevisiae"].copy()
    sdf = sdf[sdf["raw_protein_id"].astype(str).ne("")]
    return dict(zip(sdf["raw_protein_id"], sdf["canonical_gene_id"]))


def _celegans_alias_to_canonical(master_df, alias_df):
    sdf = master_df[master_df["species"] == "celegans"].copy()
    wbgene_to_canonical = dict(zip(sdf["raw_gene_id"], sdf["canonical_gene_id"]))
    alias_df = alias_df[alias_df["alias"].astype(str).ne("")].copy()
    alias_df["canonical_gene_id"] = alias_df["WBGene_id"].map(wbgene_to_canonical).fillna("")
    alias_df = alias_df[alias_df["canonical_gene_id"].astype(str).ne("")].copy()
    alias_df = alias_df.drop_duplicates(subset=["alias"], keep="first")
    return dict(zip(alias_df["alias"], alias_df["canonical_gene_id"]))


def _load_species_nodes(species):
    path = "outputs/support_graphs/{}_nodes.tsv".format(species)
    return _read_tsv(path)


def _infer_asset_rows():
    rows = []
    species_map = {
        "human": "human",
        "scerevisiae": "yeast",
        "celegans": "celegans",
    }
    assets = {
        "orthology": "Orthologs/orthologs.csv",
        "expression": "Expression/profile.csv",
        "localization": "SubLocalizations/subloc.csv",
        "prior_score": "EssentialGenes/ogee.csv",
    }
    for species, epgat_species in species_map.items():
        for group, suffix in assets.items():
            file_path = "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/{}/{}".format(epgat_species, suffix)
            exists = os.path.exists(file_path)
            source_column = ""
            usable_columns = ""
            need_remap = ""
            usable_now = "missing"
            notes = ""
            if exists:
                df = pd.read_csv(file_path, sep=None, engine="python", nrows=3).fillna("")
                columns = list(df.columns)
                source_column = "Gene" if "Gene" in columns else columns[0]
                usable_columns = ",".join(columns[: min(len(columns), 8)])
                if group == "orthology":
                    need_remap = "yes"
                    usable_now = "needs_remapping"
                    notes = "Orthology rows need species-specific raw-id to canonical_gene_id remapping before use."
                elif group == "expression":
                    need_remap = "yes"
                    usable_now = "needs_remapping"
                    notes = "Expression matrix exists but is not integrated in this round."
                elif group == "localization":
                    need_remap = "yes"
                    usable_now = "needs_remapping"
                    notes = "Localization matrix exists but is not integrated in this round."
                elif group == "prior_score":
                    need_remap = "yes"
                    usable_now = "needs_remapping"
                    notes = "OGEE label-like asset exists but is not used as a direct numeric feature in current support graph baseline."
            else:
                need_remap = "n/a"
                usable_now = "missing"
                notes = "Asset file is missing."
            rows.append(
                {
                    "species": species,
                    "feature_group": group,
                    "file_path": file_path,
                    "join_key": source_column if source_column else "unknown",
                    "usable_columns": usable_columns,
                    "need_remapping": need_remap,
                    "usable_now": usable_now,
                    "notes": notes,
                }
            )
    return pd.DataFrame(rows)


def _build_human_orthology_features(master_df):
    raw_to_canonical = _human_raw_protein_to_canonical(master_df)
    path = "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/human/Orthologs/orthologs.csv"
    df = pd.read_csv(path).fillna(0)
    feature_cols = [c for c in df.columns if c != "Gene"]
    df["canonical_gene_id"] = df["Gene"].astype(str).map(raw_to_canonical).fillna("")
    df = df[df["canonical_gene_id"].astype(str).ne("")].copy()
    if "10" in df.columns:
        support_presence = df["10"].astype(float)
    else:
        support_presence = 0.0
    out = pd.DataFrame(
        {
            "canonical_gene_id": df["canonical_gene_id"].astype(str),
            "ortholog_count": df[feature_cols].astype(float).sum(axis=1),
            "orthogroup_size": df[feature_cols].astype(float).sum(axis=1) + 1.0,
            "support_species_presence_count": support_presence,
            "conserved_across_support_species": (support_presence >= 1.0).astype(float),
            "single_copy_like": (df[feature_cols].astype(float).sum(axis=1) == 1.0).astype(float),
            "has_orthology_feature": 1.0,
        }
    )
    out = out.groupby("canonical_gene_id", as_index=False).max()
    return out, "Gene(raw_protein_id) -> canonical_gene_id via master_label_table.preliminary.raw_protein_id; support_species_presence_count uses column 10 mapped from H.sapiens-S.cerevisiae.orthoXML by EPGAT process order."


def _build_yeast_orthology_features(master_df):
    raw_to_canonical = _yeast_raw_protein_to_canonical(master_df)
    path = "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/yeast/Orthologs/orthologs.csv"
    df = pd.read_csv(path).fillna(0)
    feature_cols = [c for c in df.columns if c != "Gene"]
    df["canonical_gene_id"] = df["Gene"].astype(str).map(raw_to_canonical).fillna("")
    df = df[df["canonical_gene_id"].astype(str).ne("")].copy()
    ortho_sum = df[feature_cols].astype(float).sum(axis=1)
    out = pd.DataFrame(
        {
            "canonical_gene_id": df["canonical_gene_id"].astype(str),
            "ortholog_count": ortho_sum,
            "orthogroup_size": ortho_sum + 1.0,
            "support_species_presence_count": 0.0,
            "conserved_across_support_species": 0.0,
            "single_copy_like": (ortho_sum == 1.0).astype(float),
            "has_orthology_feature": 1.0,
        }
    )
    out = out.groupby("canonical_gene_id", as_index=False).max()
    return out, "Gene(raw_protein_id) -> canonical_gene_id via master_label_table.preliminary.raw_protein_id; current yeast orthology asset has no explicit human/celegans counterpart column, so support_species_presence_count is conservatively set to 0."


def _build_celegans_orthology_features(master_df, alias_df):
    alias_to_canonical = _celegans_alias_to_canonical(master_df, alias_df)
    path = "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/celegans/Orthologs/orthologs.csv"
    df = pd.read_csv(path).fillna(0)
    feature_cols = [c for c in df.columns if c != "Gene"]
    df["canonical_gene_id"] = df["Gene"].astype(str).map(alias_to_canonical).fillna("")
    df = df[df["canonical_gene_id"].astype(str).ne("")].copy()
    human_support = df["9606"].astype(float) if "9606" in df.columns else 0.0
    ortho_sum = df[feature_cols].astype(float).sum(axis=1)
    out = pd.DataFrame(
        {
            "canonical_gene_id": df["canonical_gene_id"].astype(str),
            "ortholog_count": ortho_sum,
            "orthogroup_size": ortho_sum + 1.0,
            "support_species_presence_count": human_support,
            "conserved_across_support_species": 0.0,
            "single_copy_like": (ortho_sum == 1.0).astype(float),
            "has_orthology_feature": 1.0,
        }
    )
    out = out.groupby("canonical_gene_id", as_index=False).max()
    return out, "Gene(symbol) -> WBGene via wb_alias_mapping.tsv -> canonical_gene_id via master_label_table.preliminary.raw_gene_id; support_species_presence_count uses taxid column 9606 when available."


def _finalize_species_matrix(species, feature_df, note):
    nodes_df = _load_species_nodes(species)
    merged = nodes_df[["species", "canonical_gene_id"]].copy().merge(
        feature_df, on="canonical_gene_id", how="left"
    )
    for col in [
        "ortholog_count",
        "orthogroup_size",
        "support_species_presence_count",
        "conserved_across_support_species",
        "single_copy_like",
        "has_orthology_feature",
    ]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0.0)
    merged["orthology_missing_mask"] = (merged["has_orthology_feature"] <= 0).astype(float)
    return merged, note


def build_all():
    outdir = "outputs/support_graph_features"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    asset_audit = _infer_asset_rows()
    asset_audit.to_csv(os.path.join(outdir, "legacy_support_feature_asset_audit.tsv"), sep="\t", index=False)

    master_df = _load_master_label_table()
    alias_df = _load_wb_alias()

    builders = {
        "human": lambda: _build_human_orthology_features(master_df),
        "scerevisiae": lambda: _build_yeast_orthology_features(master_df),
        "celegans": lambda: _build_celegans_orthology_features(master_df, alias_df),
    }

    schema_rows = []
    summary_lines = [
        "# Support Feature Integration Summary",
        "",
        "## 当前策略",
        "- 默认物种集合保持：human + scerevisiae + celegans",
        "- partial-support feature policy: missing -> 0，并保留 mask 特征",
        "- 本轮只把 orthology/conservation 真正接入支持图训练；expression / localization / prior_score 只做资产审计，不伪造接入",
        "",
        "## 每物种 orthology 集成结果",
    ]

    for species in ["human", "scerevisiae", "celegans"]:
        feature_df, note = builders[species]()
        matrix_df, note = _finalize_species_matrix(species, feature_df, note)
        matrix_path = os.path.join(outdir, "support_feature_matrix_{}.tsv".format(species))
        matrix_df.to_csv(matrix_path, sep="\t", index=False)

        coverage = float((matrix_df["has_orthology_feature"] > 0).mean()) if len(matrix_df) else 0.0
        missing_rate = 1.0 - coverage
        summary_lines.append(
            "- {}: rows={}, orthology_coverage={:.4f}, missing_rate={:.4f}".format(
                species, len(matrix_df), coverage, missing_rate
            )
        )
        summary_lines.append("  - {}".format(note))

        for feature_name in [
            "ortholog_count",
            "orthogroup_size",
            "support_species_presence_count",
            "conserved_across_support_species",
            "single_copy_like",
            "has_orthology_feature",
            "orthology_missing_mask",
        ]:
            schema_rows.append(
                {
                    "feature_name": feature_name,
                    "species_coverage": species,
                    "source_file": "/home/jiehuang/software/fungi/EPGAT/data/essential_genes/{}/Orthologs/orthologs.csv".format(
                        "human" if species == "human" else ("yeast" if species == "scerevisiae" else "celegans")
                    ),
                    "transformation": "derived from ortholog row aggregates with missing->0 and mask columns",
                    "missing_rate": missing_rate,
                }
            )

    schema_df = pd.DataFrame(schema_rows)
    schema_df.to_csv(os.path.join(outdir, "support_feature_schema.tsv"), sep="\t", index=False)
    with open(os.path.join(outdir, "support_feature_integration_summary.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary_lines))

    audit_lines = [
        "# Legacy Support Feature Asset Audit",
        "",
        "## 分类标准",
        "- directly usable: 可直接 canonical join 并立刻接入",
        "- needs remapping: 原始文件存在，但需要 species-specific namespace 转换后才能接入",
        "- unusable: 当前结构或语义不适合作为此轮 support graph 特征",
        "- missing: 文件不存在",
        "",
        "## 当前结论",
        "- orthology 是本轮优先集成对象，三物种都存在 legacy 资产，但都需要 remapping。",
        "- prior_score 当前只有 OGEE/label-like 文件，不作为 support graph 的直接 numeric prior 输入。",
        "- expression / localization 资产存在，但这轮先不接入 GraphSAGE feature matrix。",
    ]
    with open(os.path.join(outdir, "legacy_support_feature_asset_audit.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(audit_lines))
    with open("76_support_feature_asset_audit.md", "w", encoding="utf-8") as handle:
        handle.write("\n".join(audit_lines))

    print("【开始 support 特征资产审计】")
    print("当前默认物种：")
    print("  • human")
    print("  • scerevisiae")
    print("  • celegans")
    print("特征策略：")
    print("  • partial-support（允许缺失）")
    print("  • missing -> 0 + mask")
    print("orthology/conservation 特征矩阵已生成到 outputs/support_graph_features/")


if __name__ == "__main__":
    enable_utf8_stdout()
    build_all()
