import argparse
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare protocolized Fusarium label-materialization source tables"
    )
    parser.add_argument("--config", required=True, type=str)
    return parser.parse_args()


def load_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def build_old_gene_mapping(
    old_gene_list_path: str | Path,
    master_evidence_path: str | Path,
    source_mapping_path: str | Path,
    bridge_path: str | Path,
) -> pd.DataFrame:
    old = read_tsv(old_gene_list_path).drop_duplicates(subset=["Ensembl"], keep="first").copy()
    master = read_tsv(master_evidence_path)
    fg_master = master[master["species"].eq("fgraminearum")][["canonical_gene_id", "raw_gene_id"]].drop_duplicates().copy()

    raw_to_canonical = (
        fg_master.groupby("raw_gene_id")["canonical_gene_id"]
        .agg(lambda values: sorted(set(str(value).strip() for value in values if str(value).strip())))
        .to_dict()
    )

    source_mapping = read_tsv(source_mapping_path)
    source_to_canonical = (
        source_mapping.groupby("source_id")["canonical_gene_id"]
        .agg(lambda values: sorted(set(str(value).strip() for value in values if str(value).strip())))
        .to_dict()
    )

    bridge = read_tsv(bridge_path)
    xp_to_canonical = (
        bridge[bridge["bridge_status"].eq("resolved")]
        .groupby("source_protein_id")["resolved_canonical_gene_id"]
        .agg(lambda values: sorted(set(str(value).strip() for value in values if str(value).strip())))
        .to_dict()
    )

    rows: list[dict[str, str]] = []
    for _, row in old.iterrows():
        source_gene_id = str(row["Ensembl"]).strip()
        target_label = str(row["Target"]).strip()

        canonical_gene_id = ""
        mapping_status = "unresolved"
        mapping_rule = "unresolved"
        mapping_source = ""

        if source_gene_id.startswith("fgraminearum::FGRAMPH1_"):
            canonical_gene_id = source_gene_id
            mapping_status = "matched"
            mapping_rule = "exact_canonical_match"
            mapping_source = "gene_list.txt -> canonical_gene_id"
        elif source_gene_id.startswith("FGRAMPH1_"):
            canonical_gene_id = f"fgraminearum::{source_gene_id}"
            mapping_status = "matched"
            mapping_rule = "raw_gene_exact_match"
            mapping_source = "gene_list.txt -> raw_gene_id -> canonical_gene_id"
        elif source_gene_id in raw_to_canonical and len(raw_to_canonical[source_gene_id]) == 1:
            canonical_gene_id = raw_to_canonical[source_gene_id][0]
            mapping_status = "matched"
            mapping_rule = "raw_gene_exact_match"
            mapping_source = "gene_list.txt -> master_evidence raw_gene_id -> canonical_gene_id"
        elif source_gene_id in source_to_canonical and len(source_to_canonical[source_gene_id]) == 1:
            canonical_gene_id = source_to_canonical[source_gene_id][0]
            mapping_status = "matched"
            mapping_rule = "supplemental_source_match"
            mapping_source = "gene_list.txt -> protocolized source_to_canonical mapping"
        elif source_gene_id in xp_to_canonical and len(xp_to_canonical[source_gene_id]) == 1:
            canonical_gene_id = xp_to_canonical[source_gene_id][0]
            mapping_status = "matched"
            mapping_rule = "bridge_match"
            mapping_source = "gene_list.txt -> protocolized protein bridge -> canonical_gene_id"
        elif source_gene_id in raw_to_canonical and len(raw_to_canonical[source_gene_id]) > 1:
            mapping_status = "ambiguous"
            mapping_rule = "multiple_master_evidence_candidates"
            mapping_source = "gene_list.txt -> master_evidence raw_gene_id -> multiple canonical targets"
        elif source_gene_id in source_to_canonical and len(source_to_canonical[source_gene_id]) > 1:
            mapping_status = "ambiguous"
            mapping_rule = "multiple_protocolized_source_candidates"
            mapping_source = "gene_list.txt -> protocolized source_to_canonical mapping -> multiple canonical targets"
        elif source_gene_id in xp_to_canonical and len(xp_to_canonical[source_gene_id]) > 1:
            mapping_status = "ambiguous"
            mapping_rule = "multiple_protocolized_bridge_candidates"
            mapping_source = "gene_list.txt -> protocolized protein bridge -> multiple canonical targets"

        rows.append(
            {
                "source_gene_id": source_gene_id,
                "canonical_gene_id": canonical_gene_id,
                "mapping_status": mapping_status,
                "mapping_rule": mapping_rule,
                "mapping_source": mapping_source,
                "target_label": target_label,
            }
        )

    return pd.DataFrame(rows)


def build_protocolized_lethal_positive_list(
    master_evidence_path: str | Path,
    source_mapping_path: str | Path,
    bridge_path: str | Path,
) -> pd.DataFrame:
    master = read_tsv(master_evidence_path)
    source_mapping = read_tsv(source_mapping_path)
    bridge = read_tsv(bridge_path)
    fg = master[master["species"].eq("fgraminearum")].copy()
    lethal = fg[
        fg["evidence_source"].eq("phi-base_current.csv")
        & fg["evidence_term_raw"].str.lower().eq("lethal")
        & fg["supports_gold_label"].str.lower().eq("true")
        & fg["canonical_gene_id"].ne("")
    ].copy()

    # Transcript-only raw IDs were part of the historical exploratory export but were
    # not retained in the 77-gene lethal provenance table consumed by this workflow.
    lethal = lethal[~lethal["raw_gene_id"].str.match(r"^FGRAMPH1_01T\d+$", na=False)].copy()

    source_to_canonical = (
        source_mapping.groupby("source_id")["canonical_gene_id"]
        .agg(lambda values: sorted(set(str(value).strip() for value in values if str(value).strip())))
        .to_dict()
    )
    fgsg_to_canonical = (
        bridge[bridge["bridge_status"].eq("resolved")]
        .groupby("header_fgsg_id")["resolved_canonical_gene_id"]
        .agg(lambda values: sorted(set(str(value).strip() for value in values if str(value).strip())))
        .to_dict()
    )

    def resolve_final_canonical(row: pd.Series) -> str:
        canonical_gene_id = str(row["canonical_gene_id"]).strip()
        raw_gene_id = str(row["raw_gene_id"]).strip()
        if canonical_gene_id.startswith("fgraminearum::FGRAMPH1_"):
            return canonical_gene_id
        if raw_gene_id.startswith("FGRAMPH1_"):
            normalized = raw_gene_id.replace("_0G", "_01G")
            return f"fgraminearum::{normalized}"
        if raw_gene_id in source_to_canonical and len(source_to_canonical[raw_gene_id]) == 1:
            return source_to_canonical[raw_gene_id][0]
        if raw_gene_id in fgsg_to_canonical and len(fgsg_to_canonical[raw_gene_id]) == 1:
            return fgsg_to_canonical[raw_gene_id][0]
        return ""

    lethal["final_canonical_gene_id"] = lethal.apply(resolve_final_canonical, axis=1)
    lethal["final_id_status"] = lethal["final_canonical_gene_id"].map(
        lambda value: "resolved" if str(value).strip() else "unresolved"
    )
    lethal["selection_rule"] = (
        "phi_lethal_supports_gold_excluding_transcript_only_raw_ids_and_requiring_final_canonical_id"
    )
    lethal = lethal[lethal["final_canonical_gene_id"].ne("")].copy()
    lethal["canonical_gene_id"] = lethal["final_canonical_gene_id"]
    lethal = lethal.drop(columns=["final_canonical_gene_id"])
    lethal = lethal.sort_values(["canonical_gene_id", "raw_gene_id"], kind="stable").drop_duplicates(
        subset=["canonical_gene_id"], keep="first"
    )
    return lethal.reset_index(drop=True)


def build_protocolized_negative_pool(
    yeast_transfer_path: str | Path,
    bridge_path: str | Path,
    master_evidence_path: str | Path,
    source_mapping_path: str | Path,
) -> set[str]:
    yeast = read_tsv(yeast_transfer_path)
    bridge = read_tsv(bridge_path)
    evidence = read_tsv(master_evidence_path)

    yeast_bridge = yeast.merge(
        bridge[["source_protein_id", "resolved_canonical_gene_id", "bridge_status"]],
        left_on="ph1_gene_id",
        right_on="source_protein_id",
        how="left",
    )
    yeast_bridge["canonical_gene_id"] = yeast_bridge["resolved_canonical_gene_id"].fillna("")
    yeast_bridge = yeast_bridge[
        yeast_bridge["bridge_status"].eq("resolved") & yeast_bridge["canonical_gene_id"].ne("")
    ].copy()

    lethal_set = set(
        build_protocolized_lethal_positive_list(
            master_evidence_path=master_evidence_path,
            source_mapping_path=source_mapping_path,
            bridge_path=bridge_path,
        )["canonical_gene_id"].astype(str)
    )
    high_set = set(
        yeast_bridge.loc[
            yeast_bridge["weak_positive_confidence"].eq("high"), "canonical_gene_id"
        ].astype(str)
    )
    none_set = set(
        yeast_bridge.loc[
            yeast_bridge["weak_positive_confidence"].eq("none"), "canonical_gene_id"
        ].astype(str)
    )

    fg_evidence = evidence[evidence["species"].eq("fgraminearum")].copy()
    virulence_rows = fg_evidence[
        (fg_evidence["evidence_class"].astype(str) == "virulence_only")
        | fg_evidence["evidence_term_raw"].astype(str).str.contains(
            "virulence|pathogenicity", case=False, regex=True
        )
    ]
    virulence_set = set(virulence_rows["canonical_gene_id"].astype(str)) - {""}
    return set(sorted(none_set - virulence_set - lethal_set - high_set))


def build_old440_summary(
    mapping_df: pd.DataFrame,
    negative_pool: set[str],
) -> pd.DataFrame:
    positive_df = mapping_df[
        mapping_df["mapping_status"].eq("matched")
        & mapping_df["target_label"].eq("1")
        & mapping_df["canonical_gene_id"].ne("")
    ].drop_duplicates(subset=["canonical_gene_id"], keep="first")
    positive_set = set(positive_df["canonical_gene_id"].astype(str))
    overlap = positive_set & negative_pool
    final_negative_count = len(negative_pool - positive_set)
    return pd.DataFrame(
        [
            {
                "gene_list_raw_rows": int(len(mapping_df)),
                "source_gene_count_deduplicated": int(mapping_df["source_gene_id"].nunique()),
                "matched_count": int(mapping_df["mapping_status"].eq("matched").sum()),
                "unresolved_count": int(mapping_df["mapping_status"].eq("unresolved").sum()),
                "ambiguous_count": int(mapping_df["mapping_status"].eq("ambiguous").sum()),
                "final_positive_old440_count": int(len(positive_set)),
                "final_negative_old440_count": int(final_negative_count),
                "negative_overlap_removed_count": int(len(overlap)),
                "negative_pool_definition": "protocolized_newlabel_negative_pool",
            }
        ]
    )


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    paths = config["paths"]

    old_mapping = build_old_gene_mapping(
        old_gene_list_path=paths["old_gene_list"],
        master_evidence_path=paths["master_evidence_mirror"],
        source_mapping_path=Path(paths["bridge_dir"]) / "source_to_canonical_mapping.tsv",
        bridge_path=Path(paths["bridge_dir"]) / "protein_to_canonical_bridge.tsv",
    )
    lethal = build_protocolized_lethal_positive_list(
        master_evidence_path=paths["master_evidence_mirror"],
        source_mapping_path=Path(paths["bridge_dir"]) / "source_to_canonical_mapping.tsv",
        bridge_path=Path(paths["bridge_dir"]) / "protein_to_canonical_bridge.tsv",
    )
    negative_pool = build_protocolized_negative_pool(
        yeast_transfer_path=paths["yeast_transfer_table"],
        bridge_path=Path(paths["bridge_dir"]) / "protein_to_canonical_bridge.tsv",
        master_evidence_path=paths["master_evidence_mirror"],
        source_mapping_path=Path(paths["bridge_dir"]) / "source_to_canonical_mapping.tsv",
    )
    old_summary = build_old440_summary(old_mapping, negative_pool)

    Path(paths["old440_mapping_audit"]).parent.mkdir(parents=True, exist_ok=True)
    old_mapping.to_csv(paths["old440_mapping_audit"], sep="\t", index=False)
    old_summary.to_csv(paths["old440_label_summary"], sep="\t", index=False)
    lethal.to_csv(paths["lethal_positive_gene_list"], sep="\t", index=False)


if __name__ == "__main__":
    main()
