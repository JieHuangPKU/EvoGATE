import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from sklearn.model_selection import train_test_split


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Materialize Fusarium old/new label regimes into processed-data outputs")
    parser.add_argument("--config", required=True, type=str)
    return parser.parse_args()


def load_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def graph_gene_id(canonical_gene_id: str) -> str:
    value = str(canonical_gene_id).strip()
    prefix = "fgraminearum::"
    return value[len(prefix):] if value.startswith(prefix) else value


def assign_splits(labels: pd.DataFrame, seed: int, val_fraction: float, test_fraction: float, protocol_version: str) -> pd.DataFrame:
    working = labels.copy().reset_index(drop=True)
    y = working["label"].astype(int)
    train_val_idx, test_idx = train_test_split(
        working.index,
        test_size=test_fraction,
        random_state=seed,
        stratify=y,
    )
    val_relative = val_fraction / (1.0 - test_fraction)
    train_idx, val_idx = train_test_split(
        train_val_idx,
        test_size=val_relative,
        random_state=seed,
        stratify=y.loc[train_val_idx],
    )

    working["split"] = ""
    working.loc[train_idx, "split"] = "train"
    working.loc[val_idx, "split"] = "val"
    working.loc[test_idx, "split"] = "test"
    working["split_seed"] = int(seed)
    working["split_strategy"] = "stratified_fixed"
    working["test_fraction"] = float(test_fraction)
    working["val_fraction"] = float(val_fraction)
    working["split_version"] = (
        f"{protocol_version}_seed{int(seed)}_test{int(test_fraction * 100):02d}_val{int(val_fraction * 100):02d}"
    )
    return working.sort_values(["split", "canonical_gene_id"], kind="stable").reset_index(drop=True)


def ensure_no_overlap(positive_df: pd.DataFrame, negative_df: pd.DataFrame, regime: str) -> None:
    overlap = sorted(set(positive_df["canonical_gene_id"]).intersection(set(negative_df["canonical_gene_id"])))
    if overlap:
        raise ValueError(f"{regime} positive/negative overlap detected: {overlap[:10]}")


def build_newlabel(config: dict[str, Any]) -> dict[str, Any]:
    paths = config["paths"]
    bridge_dir = Path(paths["bridge_dir"])

    lethal_df = read_tsv(paths["lethal_positive_gene_list"]).drop_duplicates("canonical_gene_id", keep="first")
    yeast_df = read_tsv(paths["yeast_transfer_table"])
    evidence_df = read_tsv(paths["master_evidence_mirror"])
    bridge_df = read_tsv(bridge_dir / "protein_to_canonical_bridge.tsv")
    high_candidates_df = read_tsv(bridge_dir / "high_confidence_yeast_transfer_candidates.tsv")
    unresolved_high_df = read_tsv(bridge_dir / "unresolved_high_confidence_ids.tsv")
    bridge_summary_df = read_tsv(bridge_dir / "bridge_summary.tsv")

    lethal_set = set(lethal_df["canonical_gene_id"].astype(str))
    high_set = set(high_candidates_df["canonical_gene_id"].astype(str))
    weak_component = sorted(high_set - lethal_set)

    yeast_bridge = yeast_df.merge(
        bridge_df[["source_protein_id", "resolved_canonical_gene_id", "bridge_status", "bridge_method"]],
        left_on="ph1_gene_id",
        right_on="source_protein_id",
        how="left",
    )
    yeast_bridge["canonical_gene_id"] = yeast_bridge["resolved_canonical_gene_id"].fillna("")
    yeast_bridge = yeast_bridge[yeast_bridge["bridge_status"].eq("resolved") & yeast_bridge["canonical_gene_id"].ne("")].copy()

    none_set = set(yeast_bridge.loc[yeast_bridge["weak_positive_confidence"].eq("none"), "canonical_gene_id"].astype(str))
    fg_evidence = evidence_df[evidence_df["species"].astype(str) == "fgraminearum"].copy()
    virulence_rows = fg_evidence[
        (fg_evidence["evidence_class"].astype(str) == "virulence_only")
        | fg_evidence["evidence_term_raw"].astype(str).str.contains("virulence|pathogenicity", case=False, regex=True)
    ]
    virulence_set = set(virulence_rows["canonical_gene_id"].astype(str)) - {""}

    positive_set = sorted(lethal_set | high_set)
    negative_set = sorted(none_set - virulence_set - lethal_set - high_set)

    positive_rows = []
    transfer_support = (
        high_candidates_df.groupby("canonical_gene_id", as_index=False)
        .agg(
            supporting_xp_ids=("ph1_gene_id", lambda values: ";".join(sorted(set(map(str, values))))),
            orthogroup_ids=("orthogroup_id", lambda values: ";".join(sorted(set(map(str, values))))),
            bridge_methods=("bridge_method", lambda values: ";".join(sorted(set(map(str, values))))),
        )
        .set_index("canonical_gene_id")
        .to_dict(orient="index")
    )
    for canonical_gene_id in positive_set:
        in_lethal = canonical_gene_id in lethal_set
        in_transfer = canonical_gene_id in high_set
        if in_lethal and in_transfer:
            bucket = "lethal_plus_high_confidence_yeast_transfer"
            positive_sources = "lethal;weak_positive"
        elif in_lethal:
            bucket = "lethal_phi_supported_positive"
            positive_sources = "lethal"
        else:
            bucket = "high_confidence_yeast_transfer_positive"
            positive_sources = "weak_positive"
        support = transfer_support.get(canonical_gene_id, {})
        positive_rows.append(
            {
                "canonical_gene_id": canonical_gene_id,
                "graph_gene_id": graph_gene_id(canonical_gene_id),
                "label": 1,
                "label_text": "essential",
                "regime": "newlabel",
                "positive_sources": positive_sources,
                "construction_bucket": bucket,
                "source_manifest": str(bridge_dir / "high_confidence_yeast_transfer_candidates.tsv")
                if in_transfer and not in_lethal
                else str(paths["lethal_positive_gene_list"]),
                "support_from_lethal_positive_list": str(in_lethal).lower(),
                "support_from_protocolized_bridge": str(in_transfer).lower(),
                "supporting_xp_ids": support.get("supporting_xp_ids", ""),
                "supporting_orthogroup_ids": support.get("orthogroup_ids", ""),
                "supporting_bridge_methods": support.get("bridge_methods", ""),
                "construction_note": (
                    "Included in the rebuilt canonical newlabel positive set after protocolized bridge reconstruction."
                ),
            }
        )

    negative_rows = []
    for canonical_gene_id in negative_set:
        negative_rows.append(
            {
                "canonical_gene_id": canonical_gene_id,
                "graph_gene_id": graph_gene_id(canonical_gene_id),
                "label": 0,
                "label_text": "non-essential",
                "regime": "newlabel",
                "construction_bucket": "weak_none_after_virulence_and_positive_exclusion",
                "source_manifest": str(paths["yeast_transfer_table"]),
                "construction_note": (
                    "Included in the rebuilt canonical newlabel negative set after protocolized bridge reconstruction."
                ),
            }
        )

    positive_df = pd.DataFrame(positive_rows).sort_values("canonical_gene_id", kind="stable").reset_index(drop=True)
    negative_df = pd.DataFrame(negative_rows).sort_values("canonical_gene_id", kind="stable").reset_index(drop=True)
    labels_df = pd.concat([positive_df, negative_df], ignore_index=True)
    labels_df["source_gene_id"] = labels_df["canonical_gene_id"]
    labels_df["label_status"] = "processed_materialized"
    labels_df["label_source_project"] = "fgraminearum_newlabel"
    labels_df["label_source_file"] = labels_df["source_manifest"]

    split_df = assign_splits(
        labels=labels_df[["canonical_gene_id", "graph_gene_id", "source_gene_id", "label", "label_text"]].copy(),
        seed=int(config["runtime"]["split_seed"]),
        val_fraction=float(config["runtime"]["val_fraction"]),
        test_fraction=float(config["runtime"]["test_fraction"]),
        protocol_version=str(config["runtime"]["protocol_version"]),
    )

    audit_df = labels_df[
        [
            "canonical_gene_id",
            "graph_gene_id",
            "label",
            "label_text",
            "regime",
            "construction_bucket",
            "source_manifest",
            "label_source_project",
            "construction_note",
        ]
    ].copy()
    audit_df["is_in_lethal_positive_list"] = audit_df["canonical_gene_id"].isin(lethal_set).map(lambda value: str(value).lower())
    audit_df["is_in_protocolized_high_transfer_set"] = audit_df["canonical_gene_id"].isin(high_set).map(lambda value: str(value).lower())
    audit_df["is_in_rebuilt_negative_set"] = audit_df["canonical_gene_id"].isin(set(negative_set)).map(lambda value: str(value).lower())

    split_counts = split_df["split"].value_counts().to_dict()
    split_pos_counts = split_df[split_df["label"].astype(int) == 1]["split"].value_counts().to_dict()
    split_neg_counts = split_df[split_df["label"].astype(int) == 0]["split"].value_counts().to_dict()
    historical_positive_count = 1096
    historical_negative_count = 10270
    summary_df = pd.DataFrame(
        [
            {
                "regime": "newlabel",
                "positive_count": int(len(positive_df)),
                "negative_count": int(len(negative_df)),
                "total_count": int(len(labels_df)),
                "lethal_positive_count": int(len(lethal_set)),
                "high_confidence_yeast_transfer_component_count": int(len(weak_component)),
                "train_count": int(split_counts.get("train", 0)),
                "val_count": int(split_counts.get("val", 0)),
                "test_count": int(split_counts.get("test", 0)),
                "train_positive_count": int(split_pos_counts.get("train", 0)),
                "val_positive_count": int(split_pos_counts.get("val", 0)),
                "test_positive_count": int(split_pos_counts.get("test", 0)),
                "train_negative_count": int(split_neg_counts.get("train", 0)),
                "val_negative_count": int(split_neg_counts.get("val", 0)),
                "test_negative_count": int(split_neg_counts.get("test", 0)),
                "split_seed": int(config["runtime"]["split_seed"]),
                "split_version": split_df["split_version"].iloc[0],
                "source_lethal_positive_list": str(paths["lethal_positive_gene_list"]),
                "source_yeast_transfer_table": str(paths["yeast_transfer_table"]),
                "source_protocolized_bridge": str(bridge_dir / "protein_to_canonical_bridge.tsv"),
                "source_protocolized_high_transfer_candidates": str(bridge_dir / "high_confidence_yeast_transfer_candidates.tsv"),
                "source_master_evidence_mirror": str(paths["master_evidence_mirror"]),
                "virulence_rows_in_mirror": int(len(virulence_rows)),
                "protocolized_positive_overlap_with_lethal": int(len(high_set & lethal_set)),
                "protocolized_positive_only_yeast_transfer_count": int(len(weak_component)),
                "none_total_resolved_unique_genes": int(len(none_set)),
                "virulence_excluded_from_none": int(len(none_set & virulence_set)),
                "historical_materialized_positive_count": historical_positive_count,
                "historical_materialized_negative_count": historical_negative_count,
                "delta_vs_historical_positive_count": int(len(positive_df) - historical_positive_count),
                "delta_vs_historical_negative_count": int(len(negative_df) - historical_negative_count),
                "bridge_high_confidence_protocolized_unique_genes": int(
                    pd.to_numeric(bridge_summary_df["high_confidence_candidates_protocolized_unique_genes"], errors="coerce").fillna(0).iloc[0]
                ),
                "bridge_high_confidence_unresolved_rows": int(
                    pd.to_numeric(bridge_summary_df["high_confidence_candidates_unresolved_rows"], errors="coerce").fillna(0).iloc[0]
                ),
                "bridge_resolved_total": int(
                    pd.to_numeric(bridge_summary_df["bridge_resolved_count"], errors="coerce").fillna(0).iloc[0]
                ),
            }
        ]
    )

    summary_md = "\n".join(
        [
            "# Fusarium Newlabel Processed Label Summary",
            "",
            "## Regime Definition",
            "The `newlabel` regime is the current mainline Fusarium essential-gene protocol.",
            "It corresponds to the newer lethal plus evolution regime used for the publication-grade frozen benchmark.",
            "",
            "## Exact Source Files",
            f"- Lethal-positive provenance table: `{paths['lethal_positive_gene_list']}`",
            f"- Repo-local yeast-transfer table: `{paths['yeast_transfer_table']}`",
            f"- Protocolized protein-to-canonical bridge: `{bridge_dir / 'protein_to_canonical_bridge.tsv'}`",
            f"- Protocolized high-confidence transfer candidate table: `{bridge_dir / 'high_confidence_yeast_transfer_candidates.tsv'}`",
            f"- Mirrored evidence table used to document virulence exclusion logic: `{paths['master_evidence_mirror']}`",
            "",
            "## Construction Logic",
            "The final positive set is defined as the union of two components.",
            "First, the lethal component is taken from the preserved lethal-positive table, which records the PHI-supported direct lethal, failed-deletion, and non-viable evidence that survived canonical mapping review.",
            "Second, the transfer component is rebuilt from the repo-local PH-1 yeast-transfer table after explicit XP-to-canonical bridge reconstruction.",
            "Each PH-1 `XP_*` protein accession is mapped into the final `fgraminearum::FGRAMPH1_*` space through the protocolized bridge under `data/processed/essential_gene/fgraminearum/bridge/`.",
            "The rebuilt high-confidence transfer-positive component then retains only those PH-1 proteins whose bridge result resolves to exactly one canonical Fusarium gene.",
            "",
            "The negative set follows the same protocol definition as the mainline rebuild.",
            "It is the weak-confidence `none` pool after canonical bridging, followed by removal of genes flagged by virulence/pathogenicity evidence and removal of all genes retained in the positive set.",
            "The mirrored master evidence table is recorded here so the biological exclusion rule remains explicit in the processed artifact lineage.",
            "",
            "## How PHI Evidence Contributes",
            "PHI evidence contributes through the lethal component.",
            "Genes carrying direct lethal, failed-deletion, or non-viable mutant evidence and surviving canonical mapping review are retained in the lethal-positive provenance table and therefore enter the final positive set directly.",
            "Those genes anchor the experimentally supported essential core of the regime.",
            "",
            "## How Derived Labels Contribute",
            "Derived labels contribute through the high-confidence yeast-transfer component.",
            "Genes with strong orthology-based support from yeast essentiality transfer are first bridged from PH-1 `XP_*` protein IDs to canonical Fusarium genes, then deduplicated in canonical gene space, and finally combined with the lethal PHI-supported set.",
            "This removes the previous dependence on the historical `positive_set_P1.tsv` snapshot and makes the bridge logic explicit.",
            "",
            "## Final Counts",
            f"- positives: {len(positive_df)}",
            f"- negatives: {len(negative_df)}",
            f"- total labeled genes: {len(labels_df)}",
            f"- lethal PHI-supported positives: {len(lethal_set)}",
            f"- high-confidence yeast-transfer-supported positives retained after subtracting lethal genes: {len(weak_component)}",
            f"- unresolved high-confidence transfer rows excluded by the protocolized bridge: {len(unresolved_high_df)}",
            f"- historical snapshot positive count for comparison: {historical_positive_count}",
            f"- historical snapshot negative count for comparison: {historical_negative_count}",
            "",
            "## Split Policy",
            f"- split seed: {int(config['runtime']['split_seed'])}",
            f"- test fraction: {float(config['runtime']['test_fraction']):.2f}",
            f"- val fraction: {float(config['runtime']['val_fraction']):.2f}",
            f"- split version: `{split_df['split_version'].iloc[0]}`",
            f"- train/val/test counts: {int(split_counts.get('train', 0))}/{int(split_counts.get('val', 0))}/{int(split_counts.get('test', 0))}",
            "",
            "## Caveats and Limitations",
            "This workflow is now rebuilt from the protocolized bridge rather than from `positive_set_P1.tsv`.",
            "Because the protocolized bridge intentionally keeps only exact or auditable inferred mappings, rebuilt counts can differ from the historical snapshot when some `XP_*` rows still lack enough local evidence to resolve to a unique canonical Fusarium gene.",
            "Any unresolved transfer rows are written to the bridge audit outputs and can be revisited later when additional annotation evidence is protocolized.",
        ]
    ) + "\n"

    return {
        "regime": "newlabel",
        "positive_df": positive_df,
        "negative_df": negative_df,
        "labels_df": labels_df,
        "split_df": split_df,
        "summary_df": summary_df,
        "summary_md": summary_md,
        "audit_df": audit_df,
        "source_rows": [
            ("lethal_positive_gene_list", paths["lethal_positive_gene_list"], "PHI-backed lethal positive provenance table"),
            ("yeast_transfer_table", paths["yeast_transfer_table"], "repo-local yeast essentiality transfer source"),
            ("protocolized_bridge", bridge_dir / "protein_to_canonical_bridge.tsv", "protocolized XP-to-canonical bridge"),
            (
                "protocolized_high_transfer_candidates",
                bridge_dir / "high_confidence_yeast_transfer_candidates.tsv",
                "protocolized high-confidence yeast-transfer positive component",
            ),
            ("master_evidence_mirror", paths["master_evidence_mirror"], "mirrored evidence table documenting virulence/pathogenicity exclusions"),
        ],
    }


def build_oldlabel(config: dict[str, Any], newlabel_negative_df: pd.DataFrame) -> dict[str, Any]:
    paths = config["paths"]
    mapping_df = read_tsv(paths["old440_mapping_audit"])
    summary_input = read_tsv(paths["old440_label_summary"])

    positive_df = (
        mapping_df[
            (mapping_df["mapping_status"].astype(str) == "matched")
            & (mapping_df["target_label"].astype(str) == "1")
            & (mapping_df["canonical_gene_id"].astype(str) != "")
        ]
        .drop_duplicates("canonical_gene_id", keep="first")
        .copy()
    )
    positive_df = positive_df.sort_values("canonical_gene_id", kind="stable").reset_index(drop=True)
    positive_set = set(positive_df["canonical_gene_id"].astype(str))
    negative_base_set = set(newlabel_negative_df["canonical_gene_id"].astype(str))
    negative_set = sorted(negative_base_set - positive_set)

    negative_rows = []
    for canonical_gene_id in negative_set:
        negative_rows.append(
            {
                "canonical_gene_id": canonical_gene_id,
                "graph_gene_id": graph_gene_id(canonical_gene_id),
                "label": 0,
                "label_text": "non-essential",
                "regime": "oldlabel",
                "construction_bucket": "historical_negative_pool_minus_old440_positive_overlap",
                "source_manifest": str(paths["yeast_transfer_table"]),
                "construction_note": (
                    "Included in the preserved oldlabel negative pool after removing overlap with the old440 positive replay set."
                ),
            }
        )

    positive_rows = []
    for _, row in positive_df.iterrows():
        canonical_gene_id = str(row["canonical_gene_id"]).strip()
        positive_rows.append(
            {
                "canonical_gene_id": canonical_gene_id,
                "graph_gene_id": graph_gene_id(canonical_gene_id),
                "label": 1,
                "label_text": "essential",
                "regime": "oldlabel",
                "construction_bucket": "historical_old440_positive_replay",
                "source_manifest": str(paths["old440_mapping_audit"]),
                "source_gene_id": str(row["source_gene_id"]).strip(),
                "mapping_rule": str(row["mapping_rule"]).strip(),
                "mapping_source": str(row["mapping_source"]).strip(),
                "construction_note": (
                    "Included because the historical old440 `gene_list.txt` replay mapped this source gene to a canonical Fusarium gene with target label 1."
                ),
            }
        )

    positive_out = pd.DataFrame(positive_rows).sort_values("canonical_gene_id", kind="stable").reset_index(drop=True)
    negative_out = pd.DataFrame(negative_rows).sort_values("canonical_gene_id", kind="stable").reset_index(drop=True)
    ensure_no_overlap(positive_out, negative_out, "oldlabel")

    labels_df = pd.concat([positive_out, negative_out], ignore_index=True)
    labels_df["source_gene_id"] = labels_df["source_gene_id"].fillna(labels_df["canonical_gene_id"]) if "source_gene_id" in labels_df.columns else labels_df["canonical_gene_id"]
    labels_df["label_status"] = "processed_materialized"
    labels_df["label_source_project"] = "fgraminearum_oldlabel"
    labels_df["label_source_file"] = labels_df["source_manifest"]

    split_df = assign_splits(
        labels=labels_df[["canonical_gene_id", "graph_gene_id", "source_gene_id", "label", "label_text"]].copy(),
        seed=int(config["runtime"]["split_seed"]),
        val_fraction=float(config["runtime"]["val_fraction"]),
        test_fraction=float(config["runtime"]["test_fraction"]),
        protocol_version=str(config["runtime"]["protocol_version"]),
    )

    audit_df = labels_df[
        [
            "canonical_gene_id",
            "graph_gene_id",
            "label",
            "label_text",
            "regime",
            "construction_bucket",
            "source_manifest",
            "label_source_project",
            "construction_note",
        ]
    ].copy()
    audit_df["is_in_old440_positive_replay"] = audit_df["canonical_gene_id"].isin(positive_set).map(lambda value: str(value).lower())
    audit_df["is_in_base_negative_pool"] = audit_df["canonical_gene_id"].isin(negative_base_set).map(lambda value: str(value).lower())

    split_counts = split_df["split"].value_counts().to_dict()
    split_pos_counts = split_df[split_df["label"].astype(int) == 1]["split"].value_counts().to_dict()
    split_neg_counts = split_df[split_df["label"].astype(int) == 0]["split"].value_counts().to_dict()
    overlap_removed = int(len(negative_base_set & positive_set))
    summary_df = pd.DataFrame(
        [
            {
                "regime": "oldlabel",
                "positive_count": int(len(positive_out)),
                "negative_count": int(len(negative_out)),
                "total_count": int(len(labels_df)),
                "train_count": int(split_counts.get("train", 0)),
                "val_count": int(split_counts.get("val", 0)),
                "test_count": int(split_counts.get("test", 0)),
                "train_positive_count": int(split_pos_counts.get("train", 0)),
                "val_positive_count": int(split_pos_counts.get("val", 0)),
                "test_positive_count": int(split_pos_counts.get("test", 0)),
                "train_negative_count": int(split_neg_counts.get("train", 0)),
                "val_negative_count": int(split_neg_counts.get("val", 0)),
                "test_negative_count": int(split_neg_counts.get("test", 0)),
                "split_seed": int(config["runtime"]["split_seed"]),
                "split_version": split_df["split_version"].iloc[0],
                "source_old440_mapping_audit": str(paths["old440_mapping_audit"]),
                "source_old_gene_list": str(paths["old_gene_list"]),
                "source_negative_pool": str(paths["yeast_transfer_table"]),
                "old440_summary_positive_count": int(pd.to_numeric(summary_input["final_positive_old440_count"], errors="coerce").fillna(0).iloc[0]),
                "old440_summary_negative_count": int(pd.to_numeric(summary_input["final_negative_old440_count"], errors="coerce").fillna(0).iloc[0]),
                "old440_summary_overlap_removed": int(pd.to_numeric(summary_input["negative_overlap_removed_count"], errors="coerce").fillna(0).iloc[0]),
                "replayed_overlap_removed": overlap_removed,
            }
        ]
    )

    summary_md = "\n".join(
        [
            "# Fusarium Oldlabel Processed Label Summary",
            "",
            "## Regime Definition",
            "The `oldlabel` regime is the historical lethal plus virulence replay used for manuscript back-comparison.",
            "It is retained as an explicit legacy branch and is not the intended mainline Fusarium benchmark.",
            "",
            "## Exact Source Files",
            f"- Historical gene-list replay audit: `{paths['old440_mapping_audit']}`",
            f"- Historical old440 summary table: `{paths['old440_label_summary']}`",
            f"- Preserved historical source gene list path recorded for provenance: `{paths['old_gene_list']}`",
            f"- Base negative pool reused for the replay: `{paths['yeast_transfer_table']}` via the rebuilt newlabel none-derived negatives",
            "",
            "## Construction Logic",
            "The oldlabel positive set is reconstructed from the preserved old440 replay audit.",
            "Each retained positive gene corresponds to a source entry from the historical `gene_list.txt` whose `Target` label was 1 and whose canonical mapping was resolved successfully.",
            "The replay audit preserves the mapping rule and mapping source, allowing the processed artifact to point back to the exact historical replay lineage.",
            "",
            "The oldlabel negative set is defined as the preserved negative pool after removal of any gene that appears in the old440 positive replay set.",
            "This reproduces the old lethal plus virulence comparison regime without leaving the mainline benchmark directly dependent on the historical results directory.",
            "",
            "## Positive and Negative Definitions",
            "- Positives: genes with `target_label == 1` in the old440 replay audit after canonical mapping.",
            "- Negatives: genes from the preserved negative pool after overlap removal against the old440 positive set.",
            "",
            "## Final Counts",
            f"- positives: {len(positive_out)}",
            f"- negatives: {len(negative_out)}",
            f"- total labeled genes: {len(labels_df)}",
            f"- overlap removed from the negative pool during replay: {overlap_removed}",
            "",
            "## Split Policy",
            f"- split seed: {int(config['runtime']['split_seed'])}",
            f"- test fraction: {float(config['runtime']['test_fraction']):.2f}",
            f"- val fraction: {float(config['runtime']['val_fraction']):.2f}",
            f"- split version: `{split_df['split_version'].iloc[0]}`",
            f"- train/val/test counts: {int(split_counts.get('train', 0))}/{int(split_counts.get('val', 0))}/{int(split_counts.get('test', 0))}",
            "",
            "## Caveats and Limitations",
            "This regime is a controlled historical replay, not a newly curated essentiality protocol.",
            "Its biological meaning should be interpreted as the preserved old lethal plus virulence comparison branch rather than a current mainline recommendation.",
        ]
    ) + "\n"

    return {
        "regime": "oldlabel",
        "positive_df": positive_out,
        "negative_df": negative_out,
        "labels_df": labels_df,
        "split_df": split_df,
        "summary_df": summary_df,
        "summary_md": summary_md,
        "audit_df": audit_df,
        "source_rows": [
            ("old440_mapping_audit", paths["old440_mapping_audit"], "authoritative old440 replay audit"),
            ("old440_label_summary", paths["old440_label_summary"], "historical replay count summary"),
            ("old_gene_list", paths["old_gene_list"], "historical gene_list source preserved for provenance"),
            ("base_negative_pool", paths["yeast_transfer_table"], "rebuilt newlabel none-derived negative pool reused after overlap removal"),
        ],
    }


def write_regime_outputs(config: dict[str, Any], bundle: dict[str, Any]) -> None:
    regime = bundle["regime"]
    outdir = Path(config["paths"][f"{regime}_dir"])
    outdir.mkdir(parents=True, exist_ok=True)

    bundle["labels_df"].to_csv(outdir / "labels.tsv", sep="\t", index=False)
    bundle["positive_df"].to_csv(outdir / "positive_genes.tsv", sep="\t", index=False)
    bundle["negative_df"].to_csv(outdir / "negative_genes.tsv", sep="\t", index=False)
    bundle["split_df"].to_csv(outdir / "split.tsv", sep="\t", index=False)
    bundle["summary_df"].to_csv(outdir / "summary.tsv", sep="\t", index=False)
    (outdir / "summary.md").write_text(bundle["summary_md"], encoding="utf-8")
    bundle["audit_df"].to_csv(outdir / "label_construction_audit.tsv", sep="\t", index=False)

    source_rows = []
    for role, raw_path, description in bundle["source_rows"]:
        path = Path(raw_path)
        source_rows.append(
            {
                "regime": regime,
                "source_role": role,
                "path": str(path),
                "description": description,
                "exists": str(path.exists()).lower(),
                "sha256": sha256_file(path) if path.exists() and path.is_file() else "",
            }
        )
    pd.DataFrame(source_rows).to_csv(outdir / "source_manifest.tsv", sep="\t", index=False)

    metadata = {
        "regime": regime,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": "src.data.materialize_fgraminearum_label_regimes",
        "config_path": str(Path(config["_config_path"]).resolve()),
        "output_dir": str(outdir.resolve()),
        "split_seed": int(config["runtime"]["split_seed"]),
        "test_fraction": float(config["runtime"]["test_fraction"]),
        "val_fraction": float(config["runtime"]["val_fraction"]),
        "protocol_version": str(config["runtime"]["protocol_version"]),
        "positive_count": int(len(bundle["positive_df"])),
        "negative_count": int(len(bundle["negative_df"])),
        "total_count": int(len(bundle["labels_df"])),
    }
    (outdir / "build_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_comparison(config: dict[str, Any], old_bundle: dict[str, Any], new_bundle: dict[str, Any]) -> None:
    old_summary = old_bundle["summary_df"].iloc[0].to_dict()
    new_summary = new_bundle["summary_df"].iloc[0].to_dict()
    comparison_df = pd.DataFrame(
        [
            {
                "regime": "oldlabel",
                "positive_count": int(old_summary["positive_count"]),
                "negative_count": int(old_summary["negative_count"]),
                "total_count": int(old_summary["total_count"]),
                "train_count": int(old_summary["train_count"]),
                "val_count": int(old_summary["val_count"]),
                "test_count": int(old_summary["test_count"]),
                "definition_summary": "historical lethal plus virulence replay from old440 mapping audit",
                "positive_definition": "mapped historical gene_list Target=1 positives",
                "negative_definition": "preserved negative pool after overlap removal with old440 positives",
            },
            {
                "regime": "newlabel",
                "positive_count": int(new_summary["positive_count"]),
                "negative_count": int(new_summary["negative_count"]),
                "total_count": int(new_summary["total_count"]),
                "train_count": int(new_summary["train_count"]),
                "val_count": int(new_summary["val_count"]),
                "test_count": int(new_summary["test_count"]),
                "definition_summary": "current lethal plus evolution regime",
                "positive_definition": "lethal PHI-supported positives union high-confidence yeast-transfer-supported positives",
                "negative_definition": "weak-none pool after virulence/pathogenicity and positive exclusion",
            },
        ]
    )
    comparison_path = Path(config["paths"]["comparison_tsv"])
    comparison_path.parent.mkdir(parents=True, exist_ok=True)
    comparison_df.to_csv(comparison_path, sep="\t", index=False)

    comparison_md = "\n".join(
        [
            "# Fusarium Label Regime Comparison",
            "",
            "## Scope",
            "This file contrasts the preserved historical `oldlabel` regime with the current mainline `newlabel` regime after protocolized materialization into `data/processed`.",
            "",
            "## Count Summary",
            comparison_df.to_markdown(index=False),
            "",
            "## Key Definitional Differences",
            "- `oldlabel` is a historical replay. Its positives come from the old440 gene-list mapping audit and its negatives are the preserved negative pool after overlap removal.",
            "- `newlabel` is the current mainline regime. Its positives combine lethal PHI-supported genes with the high-confidence yeast-transfer-supported component, and its negatives follow the weak-none minus virulence/pathogenicity and positive exclusion rule preserved in the canonical materialized set.",
            "- The two regimes should remain separate because they answer different scientific questions: legacy back-comparison versus current publication-grade benchmark evaluation.",
        ]
    ) + "\n"
    Path(config["paths"]["comparison_md"]).write_text(comparison_md, encoding="utf-8")


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    config["_config_path"] = args.config

    new_bundle = build_newlabel(config)
    old_bundle = build_oldlabel(config, new_bundle["negative_df"])
    write_regime_outputs(config, old_bundle)
    write_regime_outputs(config, new_bundle)
    write_comparison(config, old_bundle, new_bundle)


if __name__ == "__main__":
    main()
