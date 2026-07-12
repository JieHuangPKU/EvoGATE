#!/usr/bin/env python3
"""
脚本名称: build_figure1c_label_source_sankey_data.py
日期: 2026-05-06
作者: OpenAI/Codex + Jie Huang workflow support
功能描述: 整理 Fusarium graminearum new label 正标签来源与 yeast-to-Fusarium
          transfer path，输出 Figure 1C ggsankey 绘图数据、summary box 和审计表。
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "results" / "Figure1C"

NEWLABEL_DIR = REPO_ROOT / "data/processed/essential_gene/fgraminearum/newlabel"
OLDLABEL_DIR = REPO_ROOT / "data/processed/essential_gene/fgraminearum/oldlabel"
BRIDGE_DIR = REPO_ROOT / "data/processed/essential_gene/fgraminearum/bridge"
DERIVED_DIR = REPO_ROOT / "data/derived_labels"
INTERIM_DIR = REPO_ROOT / "data/interim/protocol_refactor/fgraminearum_label_materialization"

NEW_LABELS_TSV = NEWLABEL_DIR / "labels.tsv"
NEW_POSITIVE_TSV = NEWLABEL_DIR / "positive_genes.tsv"
NEW_LABEL_AUDIT_TSV = NEWLABEL_DIR / "label_construction_audit.tsv"
NEW_SUMMARY_TSV = NEWLABEL_DIR / "summary.tsv"
OLD_LABELS_TSV = OLDLABEL_DIR / "labels.tsv"
OLD_SUMMARY_TSV = OLDLABEL_DIR / "summary.tsv"

HIGH_CONFIDENCE_TSV = BRIDGE_DIR / "high_confidence_yeast_transfer_candidates.tsv"
UNRESOLVED_HIGH_CONFIDENCE_TSV = BRIDGE_DIR / "unresolved_high_confidence_ids.tsv"
PROTEIN_BRIDGE_TSV = BRIDGE_DIR / "protein_to_canonical_bridge.tsv"
BRIDGE_SUMMARY_TSV = BRIDGE_DIR / "bridge_summary.tsv"
YEAST_TRANSFER_TSV = DERIVED_DIR / "ph1_yeast_essential_ortholog_labels.tsv"
LETHAL_POSITIVE_TSV = INTERIM_DIR / "lethal_positive_gene_list.tsv"

SOURCE_LABELS = {
    "scer_only": "S. cerevisiae only",
    "spom_only": "S. pombe only",
    "both": "Shared by both yeasts",
    "phi": "PHI-base essential/lethal",
    "other": "Other / remaining evidence",
}

SOURCE_ORDER = [
    "S. cerevisiae only",
    "S. pombe only",
    "Shared by both yeasts",
    "PHI-base essential/lethal",
    "Other / remaining evidence",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Figure 1C source-resolved transfer Sankey data.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Figure 1C output directory.")
    return parser.parse_args()


def log(message: str) -> None:
    print(f"[Figure1C data] {message}")


def require_files(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input files:\n- " + "\n- ".join(missing))


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def clean_count(value: object) -> int:
    text = str(value).strip()
    return int(float(text)) if text else 0


def source_label(raw_source: str) -> str:
    return SOURCE_LABELS.get(raw_source, SOURCE_LABELS["other"])


def compact_stage_node(source: str, stage: str, count: int | None = None) -> str:
    suffix = f" (n={count})" if count is not None else ""
    if stage == "orthogroup":
        if source == "phi":
            return "Direct PHI evidence\n(no yeast orthogroup)" + suffix
        if source == "other":
            return "Other evidence group" + suffix
        return f"{source_label(source)}\nsupported orthogroups" + suffix
    if stage == "gene":
        if source == "phi":
            return "PHI-mapped Fusarium genes" + suffix
        if source == "other":
            return "Other mapped Fusarium genes" + suffix
        return f"{source_label(source)}\nmapped Fusarium genes" + suffix
    raise ValueError(f"Unexpected compact stage: {stage}")


def load_inputs() -> dict[str, pd.DataFrame]:
    # 读取 canonical final label 数据
    return {
        "new_labels": read_tsv(NEW_LABELS_TSV),
        "new_positive": read_tsv(NEW_POSITIVE_TSV),
        "new_audit": read_tsv(NEW_LABEL_AUDIT_TSV),
        "new_summary": read_tsv(NEW_SUMMARY_TSV),
        "old_labels": read_tsv(OLD_LABELS_TSV),
        "old_summary": read_tsv(OLD_SUMMARY_TSV),
        "high_conf": read_tsv(HIGH_CONFIDENCE_TSV),
        "unresolved_high": read_tsv(UNRESOLVED_HIGH_CONFIDENCE_TSV),
        "protein_bridge": read_tsv(PROTEIN_BRIDGE_TSV),
        "yeast_transfer": read_tsv(YEAST_TRANSFER_TSV),
        "lethal": read_tsv(LETHAL_POSITIVE_TSV),
        "bridge_summary": read_tsv(BRIDGE_SUMMARY_TSV),
    }


def build_positive_source_audit(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    # 整理 yeast transfer 来源归属
    positives = tables["new_positive"].copy()
    high_conf = tables["high_conf"].copy()
    lethal = tables["lethal"].copy()
    old_positive_set = set(tables["old_labels"].loc[tables["old_labels"]["label"].eq("1"), "canonical_gene_id"])
    phi_set = set(lethal["canonical_gene_id"])

    high_conf_lookup = (
        high_conf.sort_values(["canonical_gene_id", "ph1_gene_id"], kind="stable")
        .groupby("canonical_gene_id", as_index=False)
        .agg(
            yeast_support_class=("yeast_essential_support_class", "first"),
            supporting_orthogroup_ids=("orthogroup_id", lambda values: ";".join(sorted(set(values)))),
            supporting_xp_ids=("ph1_gene_id", lambda values: ";".join(sorted(set(values)))),
            scer_essential_gene_ids=("scer_essential_gene_ids", lambda values: ";".join(sorted({v for v in values if v}))),
            spom_essential_gene_ids=("spom_essential_gene_ids", lambda values: ";".join(sorted({v for v in values if v}))),
            bridge_methods=("bridge_method", lambda values: ";".join(sorted({v for v in values if v}))),
        )
    )

    audit = positives.merge(high_conf_lookup, on="canonical_gene_id", how="left", suffixes=("", "_bridge")).fillna("")
    for column in ["supporting_orthogroup_ids", "supporting_xp_ids", "bridge_methods"]:
        bridge_column = f"{column}_bridge"
        if bridge_column in audit.columns:
            audit[column] = audit[bridge_column].where(audit[bridge_column].ne(""), audit.get(column, ""))
    audit["has_yeast_transfer_support"] = audit["support_from_protocolized_bridge"].map(yes)
    audit["has_phi_essential_lethal_support"] = audit["canonical_gene_id"].isin(phi_set)
    audit["is_old_positive"] = audit["canonical_gene_id"].isin(old_positive_set)

    # 合并 PHI-base 证据并采用主来源唯一归属，避免 Sankey double count
    def assign_primary(row: pd.Series) -> str:
        if row["has_yeast_transfer_support"]:
            raw = row["yeast_support_class"] or "other"
            return raw if raw in {"scer_only", "spom_only", "both"} else "other"
        if row["has_phi_essential_lethal_support"]:
            return "phi"
        return "other"

    audit["primary_source_key"] = audit.apply(assign_primary, axis=1)
    audit["support_source"] = audit["primary_source_key"].map(source_label)
    audit["orthogroup_stage_detail"] = audit.apply(
        lambda row: row["supporting_orthogroup_ids"]
        if row["supporting_orthogroup_ids"]
        else ("Direct PHI evidence (no yeast orthogroup)" if row["primary_source_key"] == "phi" else "Other evidence"),
        axis=1,
    )
    audit["mapped_fusarium_gene"] = audit["canonical_gene_id"]
    audit["final_label_class"] = "Essential (positive)"
    audit["supporting_evidence_flags"] = audit.apply(
        lambda row: ";".join(
            flag
            for flag, present in [
                ("yeast_transfer", row["has_yeast_transfer_support"]),
                ("phi_base_essential_lethal", row["has_phi_essential_lethal_support"]),
                ("old_positive", row["is_old_positive"]),
            ]
            if present
        ),
        axis=1,
    )

    audit_cols = [
        "canonical_gene_id",
        "graph_gene_id",
        "label",
        "construction_bucket",
        "support_source",
        "primary_source_key",
        "yeast_support_class",
        "has_yeast_transfer_support",
        "has_phi_essential_lethal_support",
        "is_old_positive",
        "supporting_orthogroup_ids",
        "orthogroup_stage_detail",
        "supporting_xp_ids",
        "scer_essential_gene_ids",
        "spom_essential_gene_ids",
        "bridge_methods",
        "mapped_fusarium_gene",
        "final_label_class",
        "supporting_evidence_flags",
        "source_manifest",
    ]
    return audit[audit_cols].sort_values(["support_source", "canonical_gene_id"], kind="stable")


def build_flow_paths(source_audit: pd.DataFrame, tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    # 构造 ggsankey 长表所需的 source-resolved transfer path
    positive_paths = []
    for source_key, group in source_audit.groupby("primary_source_key", sort=False):
        count = len(group)
        if count == 0:
            continue
        og_count = group["supporting_orthogroup_ids"].replace("", pd.NA).dropna().nunique()
        gene_count = group["canonical_gene_id"].nunique()
        positive_paths.append(
            {
                "path_id": f"{source_key}_positive",
                "support_source": source_label(source_key),
                "source_key": source_key,
                "orthogroup_stage": compact_stage_node(source_key, "orthogroup", og_count if source_key != "phi" else count),
                "mapped_gene_stage": compact_stage_node(source_key, "gene", gene_count),
                "final_label_class": "Essential (positive)",
                "value": count,
                "path_scope": "final_positive_primary_source",
            }
        )

    # 补充 high-confidence yeast transfer 中未解析/未保留的路径，提供 Non-essential / excluded 终点
    unresolved = tables["unresolved_high"].copy()
    excluded_paths = []
    if not unresolved.empty:
        unresolved["source_key"] = unresolved["yeast_essential_support_class"].where(
            unresolved["yeast_essential_support_class"].isin(["scer_only", "spom_only", "both"]), "other"
        )
        for source_key, group in unresolved.groupby("source_key", sort=False):
            count = len(group)
            og_count = group["orthogroup_id"].nunique()
            excluded_paths.append(
                {
                    "path_id": f"{source_key}_excluded_unresolved",
                    "support_source": source_label(source_key),
                    "source_key": source_key,
                    "orthogroup_stage": compact_stage_node(source_key, "orthogroup", og_count),
                    "mapped_gene_stage": "Unresolved bridge IDs\n(not retained)",
                    "final_label_class": "Non-essential / excluded",
                    "value": count,
                    "path_scope": "high_confidence_unresolved_excluded",
                }
            )

    paths = pd.DataFrame(positive_paths + excluded_paths)
    source_rank = {label: idx for idx, label in enumerate(SOURCE_ORDER)}
    paths["source_rank"] = paths["support_source"].map(source_rank).fillna(99).astype(int)
    return paths.sort_values(["source_rank", "final_label_class", "path_id"], kind="stable").reset_index(drop=True)


def paths_to_sankey_long(paths: pd.DataFrame) -> pd.DataFrame:
    rows = []
    transitions = [
        ("Support source", "support_source", "Supported orthogroups", "orthogroup_stage"),
        ("Supported orthogroups", "orthogroup_stage", "Mapped Fusarium genes", "mapped_gene_stage"),
        ("Mapped Fusarium genes", "mapped_gene_stage", "Final label class", "final_label_class"),
    ]
    for _, path in paths.iterrows():
        for stage_pair, (x, node_col, next_x, next_node_col) in enumerate(transitions, start=1):
            rows.append(
                {
                    "path_id": path["path_id"],
                    "stage_pair": stage_pair,
                    "x": x,
                    "node": path[node_col],
                    "next_x": next_x,
                    "next_node": path[next_node_col],
                    "value": int(path["value"]),
                    "support_source": path["support_source"],
                    "source_key": path["source_key"],
                    "final_label_class": path["final_label_class"],
                    "path_scope": path["path_scope"],
                }
            )
    return pd.DataFrame(rows)


def build_stage_counts(paths: pd.DataFrame) -> pd.DataFrame:
    # 输出每层节点计数，供图注和 README 复核
    records = []
    layer_cols = [
        ("Support source", "support_source"),
        ("Supported orthogroups", "orthogroup_stage"),
        ("Mapped Fusarium genes", "mapped_gene_stage"),
        ("Final label class", "final_label_class"),
    ]
    for layer, col in layer_cols:
        grouped = paths.groupby(col, as_index=False)["value"].sum().rename(columns={col: "node", "value": "count"})
        grouped["layer"] = layer
        records.append(grouped[["layer", "node", "count"]])
    return pd.concat(records, ignore_index=True).sort_values(["layer", "node"], kind="stable")


def build_summary_box(source_audit: pd.DataFrame, tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    # 输出 summary box 统计
    new_summary = tables["new_summary"].iloc[0]
    old_summary = tables["old_summary"].iloc[0]
    high_conf = tables["high_conf"]
    phi_supported = int(source_audit["has_phi_essential_lethal_support"].sum())
    yeast_supported = int(source_audit["has_yeast_transfer_support"].sum())
    summary_rows = [
        ("old positives", clean_count(old_summary["positive_count"]), rel(OLD_SUMMARY_TSV), "oldlabel positive_count"),
        ("new positives", clean_count(new_summary["positive_count"]), rel(NEW_SUMMARY_TSV), "newlabel positive_count"),
        (
            "high-confidence transferred positives",
            int(high_conf["canonical_gene_id"].nunique()),
            rel(HIGH_CONFIDENCE_TSV),
            "unique canonical_gene_id in protocolized high-confidence transfer candidates",
        ),
        ("old total", clean_count(old_summary["total_count"]), rel(OLD_SUMMARY_TSV), "oldlabel total_count"),
        ("new total", clean_count(new_summary["total_count"]), rel(NEW_SUMMARY_TSV), "newlabel total_count"),
        (
            "PHI-supported positives",
            phi_supported,
            rel(LETHAL_POSITIVE_TSV),
            "final positives with PHI-base essential/lethal support",
        ),
        (
            "yeast-transfer-supported positives",
            yeast_supported,
            rel(NEW_POSITIVE_TSV),
            "final positives with support_from_protocolized_bridge=true",
        ),
    ]
    return pd.DataFrame(summary_rows, columns=["metric", "count", "source_file", "definition"])


def build_readme(paths: pd.DataFrame, stage_counts: pd.DataFrame, summary: pd.DataFrame) -> str:
    # 写清楚 Figure 1C 统计口径
    source_counts = (
        paths.loc[paths["final_label_class"].eq("Essential (positive)")]
        .groupby("support_source")["value"]
        .sum()
        .sort_index()
    )
    source_lines = "\n".join(f"- {source}: {int(count)}" for source, count in source_counts.items())
    summary_lines = "\n".join(f"- {row.metric}: {int(row.count)}" for row in summary.itertuples())
    stage_lines = "\n".join(f"- {row.layer} | {row.node}: {int(row.count)}" for row in stage_counts.itertuples())
    used_files = [
        NEW_LABELS_TSV,
        NEW_POSITIVE_TSV,
        NEW_LABEL_AUDIT_TSV,
        NEW_SUMMARY_TSV,
        OLD_LABELS_TSV,
        OLD_SUMMARY_TSV,
        HIGH_CONFIDENCE_TSV,
        UNRESOLVED_HIGH_CONFIDENCE_TSV,
        PROTEIN_BRIDGE_TSV,
        BRIDGE_SUMMARY_TSV,
        YEAST_TRANSFER_TSV,
        LETHAL_POSITIVE_TSV,
    ]
    file_lines = "\n".join(f"- `{rel(path)}`" for path in used_files)
    return f"""# Figure 1C Source-Resolved Transfer Sankey

Figure 1C is a source-resolved transfer Sankey for the `fgraminearum` new label positive set. It is designed to answer both source composition and transfer path questions: where the final positives came from, and how yeast essential evidence passed through supported orthogroups and canonical Fusarium mapping into final retained positives.

## Actual Input Files
{file_lines}

## Canonical Final Label Source
`{rel(NEW_POSITIVE_TSV)}` is used as the canonical final positive source because it is the materialized newlabel positive subset with final retained canonical IDs, construction buckets, supporting XP IDs, orthogroup IDs, and bridge methods. `{rel(NEW_LABELS_TSV)}` is used for the canonical newlabel total, and `{rel(NEW_LABEL_AUDIT_TSV)}` is retained as the construction audit table that records the final label-construction decision fields. The summary box totals come from `{rel(NEW_SUMMARY_TSV)}` and `{rel(OLD_SUMMARY_TSV)}`.

## Source Assignment
The Sankey uses a primary-source unique assignment, so path counts close and no positive gene is silently double counted.

- `S. cerevisiae only`: final positive has `support_from_protocolized_bridge=true` and its high-confidence bridge row has `yeast_essential_support_class=scer_only`.
- `S. pombe only`: final positive has `support_from_protocolized_bridge=true` and `yeast_essential_support_class=spom_only`.
- `Shared by both yeasts`: final positive has `support_from_protocolized_bridge=true` and `yeast_essential_support_class=both`.
- `PHI-base essential/lethal`: final positive has no protocolized yeast-transfer primary support but is present in `{rel(LETHAL_POSITIVE_TSV)}`.
- Yeast + PHI overlap genes are assigned to the yeast primary source, while `figure1c_source_audit.tsv` records `has_phi_essential_lethal_support=true`.

## Transfer Path Definition
- Layer 1 `Support source`: primary source category above.
- Layer 2 `Supported orthogroups`: yeast sources use the real `orthogroup_id` support summarized into source-specific orthogroup pools for plotting; PHI-only positives use `Direct PHI evidence (no yeast orthogroup)`.
- Layer 3 `Mapped Fusarium genes`: canonical `fgraminearum::FGRAMPH1_*` genes retained in `{rel(NEW_POSITIVE_TSV)}`; unresolved high-confidence yeast rows are displayed as `Unresolved bridge IDs (not retained)` and terminate at `Non-essential / excluded`.
- Layer 4 `Final label class`: `Essential (positive)` for retained positives and `Non-essential / excluded` for unresolved high-confidence transfer candidates that did not enter the final positive set.

## Positive Source Composition
{source_lines}

## Summary Box Values
{summary_lines}

## Stage Counts
{stage_lines}

## Outputs
- `figure1c_sankey_long.tsv`: ggsankey-ready transition table with `x`, `node`, `next_x`, `next_node`, `value`, and source metadata.
- `figure1c_stage_counts.tsv`: node counts for each plotted layer.
- `figure1c_summary_box.tsv`: right-side summary box statistics.
- `figure1c_source_audit.tsv`: one row per final positive gene, preserving yeast support class, PHI support, old positive overlap, orthogroup IDs, XP IDs, and bridge method.
"""


def main() -> None:
    args = parse_args()
    require_files(
        [
            NEW_LABELS_TSV,
            NEW_POSITIVE_TSV,
            NEW_LABEL_AUDIT_TSV,
            NEW_SUMMARY_TSV,
            OLD_LABELS_TSV,
            OLD_SUMMARY_TSV,
            HIGH_CONFIDENCE_TSV,
            UNRESOLVED_HIGH_CONFIDENCE_TSV,
            PROTEIN_BRIDGE_TSV,
            BRIDGE_SUMMARY_TSV,
            YEAST_TRANSFER_TSV,
            LETHAL_POSITIVE_TSV,
        ]
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    log("读取 canonical final label、bridge、yeast transfer 与 PHI-base 数据")
    tables = load_inputs()

    log("整理 final positive gene 的主来源归属与 supporting evidence 审计字段")
    source_audit = build_positive_source_audit(tables)

    log("构造 source-resolved transfer Sankey 聚合路径")
    paths = build_flow_paths(source_audit, tables)
    sankey_long = paths_to_sankey_long(paths)
    stage_counts = build_stage_counts(paths)

    log("输出 summary box 统计和 README")
    summary_box = build_summary_box(source_audit, tables)
    readme = build_readme(paths, stage_counts, summary_box)

    sankey_path = args.output_dir / "figure1c_sankey_long.tsv"
    stage_path = args.output_dir / "figure1c_stage_counts.tsv"
    summary_path = args.output_dir / "figure1c_summary_box.tsv"
    audit_path = args.output_dir / "figure1c_source_audit.tsv"
    readme_path = args.output_dir / "README.md"

    sankey_long.to_csv(sankey_path, sep="\t", index=False)
    stage_counts.to_csv(stage_path, sep="\t", index=False)
    summary_box.to_csv(summary_path, sep="\t", index=False)
    source_audit.to_csv(audit_path, sep="\t", index=False)
    readme_path.write_text(readme, encoding="utf-8")

    log(f"Figure 1C Sankey data written to {rel(sankey_path)}")
    log(f"Figure 1C stage counts written to {rel(stage_path)}")
    log(f"Figure 1C summary box written to {rel(summary_path)}")
    log(f"Figure 1C source audit written to {rel(audit_path)}")
    log(f"Figure 1C README written to {rel(readme_path)}")
    log("Done.")


if __name__ == "__main__":
    main()
