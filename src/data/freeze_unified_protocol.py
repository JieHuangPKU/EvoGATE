import argparse
from pathlib import Path
from typing import Any, Dict, List, Union

import pandas as pd
import yaml
from sklearn.model_selection import train_test_split


STANDARD_LABEL_COLUMNS = {
    "gene_id",
    "gene_id_source",
    "gene_symbol",
    "species",
    "label_binary",
    "label_text",
    "label_source_project",
    "label_source_file",
    "evidence_note",
    "included_in_final",
    "exclusion_reason",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Freeze the publication-grade ProGATE_v2 protocol manifests")
    parser.add_argument("--config", required=True, type=str)
    return parser.parse_args()


def load_config(path: Union[str, Path]) -> Dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _read_tsv(path: Union[str, Path]) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Required TSV file not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def resolve_source_path(config: Dict[str, Any], source_key: str) -> Path:
    try:
        return Path(config["label_sources"][source_key])
    except KeyError as exc:
        raise KeyError(f"Undefined label source '{source_key}' in configs/frozen_protocol.yaml") from exc


def canonical_to_graph_gene_id(canonical_gene_id: str, species: str) -> str:
    value = str(canonical_gene_id).strip()
    prefix = f"{species}::"
    return value[len(prefix):] if value.startswith(prefix) else value


def load_standard_processed_labels(config: Dict[str, Any], protocol_name: str, protocol_cfg: Dict[str, Any]) -> pd.DataFrame:
    source_path = resolve_source_path(config, protocol_cfg["label_source_key"])
    labels = _read_tsv(source_path)
    missing = sorted(STANDARD_LABEL_COLUMNS.difference(labels.columns))
    if missing:
        raise ValueError(f"Standard processed labels missing columns {missing}: {source_path}")

    labels = labels[labels["included_in_final"].astype(str).str.lower() == "true"].copy()
    labels["gene_id"] = labels["gene_id"].astype(str).str.strip()
    labels = labels[labels["gene_id"].ne("")].drop_duplicates(subset=["gene_id"], keep="first").reset_index(drop=True)
    labels["label"] = pd.to_numeric(labels["label_binary"], errors="raise").astype(int)

    out = pd.DataFrame(
        {
            "protocol_name": protocol_name,
            "species": protocol_cfg["species"],
            "regime": protocol_cfg["regime"],
            "data_key": protocol_cfg["data_key"],
            "canonical_gene_id": labels["gene_id"],
            "graph_gene_id": labels["gene_id"],
            "source_gene_id": labels["gene_id"],
            "label": labels["label"],
            "label_text": labels["label_text"],
            "label_status": "frozen",
            "label_source_project": labels["label_source_project"],
            "label_source_file": labels["label_source_file"],
            "source_manifest": str(source_path),
            "positive_definition": labels["label_text"].map(
                lambda value: "essential" if str(value).strip().lower() == "essential" else "non-essential"
            ),
            "notes": labels["evidence_note"],
            "is_mainline_protocol": bool(protocol_cfg["is_mainline"]),
            "protocol_version": config["runtime"]["protocol_version"],
        }
    )
    return out.sort_values(["label", "canonical_gene_id"], ascending=[False, True], kind="stable").reset_index(drop=True)


def load_pair_labels(config: Dict[str, Any], protocol_name: str, protocol_cfg: Dict[str, Any]) -> pd.DataFrame:
    positive_path = resolve_source_path(config, protocol_cfg["positive_source_key"])
    negative_path = resolve_source_path(config, protocol_cfg["negative_source_key"])
    positive_df = _read_tsv(positive_path)
    negative_df = _read_tsv(negative_path)
    if "canonical_gene_id" not in positive_df.columns or "canonical_gene_id" not in negative_df.columns:
        raise ValueError(
            f"Positive/negative label tables must both contain canonical_gene_id: {positive_path}, {negative_path}"
        )

    positive_df["canonical_gene_id"] = positive_df["canonical_gene_id"].astype(str).str.strip()
    negative_df["canonical_gene_id"] = negative_df["canonical_gene_id"].astype(str).str.strip()
    positive_df = positive_df[positive_df["canonical_gene_id"].ne("")].drop_duplicates("canonical_gene_id", keep="first")
    negative_df = negative_df[negative_df["canonical_gene_id"].ne("")].drop_duplicates("canonical_gene_id", keep="first")

    overlap = sorted(set(positive_df["canonical_gene_id"]).intersection(set(negative_df["canonical_gene_id"])))
    if overlap:
        raise ValueError(
            f"Overlapping positive/negative genes detected for protocol '{protocol_name}'; first examples: {overlap[:10]}"
        )

    if protocol_cfg["regime"] == "oldlabel":
        positive_definition = "historical old-label positive set reconstructed from gene_list old440 mapping"
        negative_definition = "historical old-label negative set after overlap removal"
        label_source_project = "historical_fusarium_oldlabel"
    else:
        positive_definition = "phase2b_new_label positive set: lethal union high-confidence yeast-transfer positives (P1)"
        negative_definition = "phase2b_new_label negative set: weak-confidence none minus virulence and positive genes"
        label_source_project = "phase2b_new_label"

    positive_rows = pd.DataFrame(
        {
            "protocol_name": protocol_name,
            "species": protocol_cfg["species"],
            "regime": protocol_cfg["regime"],
            "data_key": protocol_cfg["data_key"],
            "canonical_gene_id": positive_df["canonical_gene_id"],
            "graph_gene_id": positive_df["canonical_gene_id"].map(
                lambda value: canonical_to_graph_gene_id(value, protocol_cfg["species"])
            ),
            "source_gene_id": positive_df.get("source_gene_id", positive_df["canonical_gene_id"]),
            "label": 1,
            "label_text": "essential",
            "label_status": "frozen",
            "label_source_project": label_source_project,
            "label_source_file": str(positive_path),
            "source_manifest": str(positive_path),
            "positive_definition": positive_definition,
            "notes": positive_df.get("positive_sources", ""),
            "is_mainline_protocol": bool(protocol_cfg["is_mainline"]),
            "protocol_version": config["runtime"]["protocol_version"],
        }
    )
    negative_rows = pd.DataFrame(
        {
            "protocol_name": protocol_name,
            "species": protocol_cfg["species"],
            "regime": protocol_cfg["regime"],
            "data_key": protocol_cfg["data_key"],
            "canonical_gene_id": negative_df["canonical_gene_id"],
            "graph_gene_id": negative_df["canonical_gene_id"].map(
                lambda value: canonical_to_graph_gene_id(value, protocol_cfg["species"])
            ),
            "source_gene_id": negative_df.get("source_gene_id", negative_df["canonical_gene_id"]),
            "label": 0,
            "label_text": "non-essential",
            "label_status": "frozen",
            "label_source_project": label_source_project,
            "label_source_file": str(negative_path),
            "source_manifest": str(negative_path),
            "positive_definition": negative_definition,
            "notes": negative_df.get("positive_sources", ""),
            "is_mainline_protocol": bool(protocol_cfg["is_mainline"]),
            "protocol_version": config["runtime"]["protocol_version"],
        }
    )
    out = pd.concat([positive_rows, negative_rows], ignore_index=True)
    return out.sort_values(["label", "canonical_gene_id"], ascending=[False, True], kind="stable").reset_index(drop=True)


def build_label_manifest(config: Dict[str, Any], protocol_name: str) -> pd.DataFrame:
    protocol_cfg = config["protocols"][protocol_name]
    source_type = protocol_cfg["label_source_type"]
    if source_type == "standard_processed":
        return load_standard_processed_labels(config, protocol_name, protocol_cfg)
    if source_type == "positive_negative_pair":
        return load_pair_labels(config, protocol_name, protocol_cfg)
    raise ValueError(f"Unsupported label_source_type '{source_type}' for protocol '{protocol_name}'")


def assign_splits(labels: pd.DataFrame, seed: int, val_fraction: float, test_fraction: float) -> pd.DataFrame:
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
        f"{working['protocol_version'].iloc[0]}_seed{int(seed)}_test{int(test_fraction * 100):02d}_val{int(val_fraction * 100):02d}"
    )
    return working.sort_values(["split", "canonical_gene_id"], kind="stable").reset_index(drop=True)


def label_output_path(config: Dict[str, Any], protocol_name: str) -> Path:
    return Path(config["paths"]["labels_dir"]) / str(config["protocols"][protocol_name]["label_output"])


def split_output_path(config: Dict[str, Any], protocol_name: str) -> Path:
    return Path(config["paths"]["splits_dir"]) / str(config["protocols"][protocol_name]["split_output"])


def summarize_label_counts(labels: pd.DataFrame) -> Dict[str, int]:
    return {
        "positives": int((labels["label"].astype(int) == 1).sum()),
        "negatives": int((labels["label"].astype(int) == 0).sum()),
        "total": int(len(labels)),
    }


def summarize_split_counts(split_df: pd.DataFrame) -> Dict[str, int]:
    labels = split_df["label"].astype(int)
    splits = split_df["split"].astype(str)
    return {
        "train": int((splits == "train").sum()),
        "val": int((splits == "val").sum()),
        "test": int((splits == "test").sum()),
        "train_pos": int(((splits == "train") & (labels == 1)).sum()),
        "val_pos": int(((splits == "val") & (labels == 1)).sum()),
        "test_pos": int(((splits == "test") & (labels == 1)).sum()),
        "train_neg": int(((splits == "train") & (labels == 0)).sum()),
        "val_neg": int(((splits == "val") & (labels == 0)).sum()),
        "test_neg": int(((splits == "test") & (labels == 0)).sum()),
    }


def write_label_protocol_summary(config: Dict[str, Any], label_rows: List[Dict[str, Any]]) -> None:
    labels_dir = Path(config["paths"]["labels_dir"])
    lines = [
        "# Frozen Label Protocol Summary",
        "",
        f"- protocol_version: `{config['runtime']['protocol_version']}`",
        "- benchmark species set: human, celegans, scerevisiae, dmelanogaster, fgraminearum",
        "- Fusarium mainline regime: `fgraminearum_newlabel`",
        "- Fusarium legacy replay regime: `fgraminearum_oldlabel`",
        "- Deprecated mainline inputs: `broad79`, `strict29`, `conflict8`, `fgraminearum_gold_positive*`, implicit `old440` defaults",
        "",
        "| protocol | species | regime | positives | negatives | total | source_manifest | output | mainline |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in label_rows:
        lines.append(
            f"| {row['protocol_name']} | {row['species']} | {row['regime']} | {row['positives']} | {row['negatives']} | "
            f"{row['total']} | `{row['source_manifest']}` | `{row['output_file']}` | {str(bool(row['is_mainline'])).lower()} |"
        )
    (labels_dir / "label_protocol_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_fusarium_label_summary(config: Dict[str, Any], label_rows: List[Dict[str, Any]]) -> None:
    labels_dir = Path(config["paths"]["labels_dir"])
    old_row = next(row for row in label_rows if row["protocol_name"] == "fgraminearum_oldlabel")
    new_row = next(row for row in label_rows if row["protocol_name"] == "fgraminearum_newlabel")
    lines = [
        "# Fusarium Label Protocol Summary",
        "",
        "## Mainline Decision",
        "- mainline regime: `fgraminearum_newlabel`",
        "- legacy comparison regime: `fgraminearum_oldlabel`",
        "- old440 is retained only as explicit historical replay, never as an implicit default",
        "",
        "## Old Regime",
        f"- protocol: `{old_row['protocol_name']}`",
        "- role: legacy replay / controlled historical comparison",
        f"- positive source: `{config['label_sources']['fgraminearum_oldlabel_positive']}`",
        f"- negative source: `{config['label_sources']['fgraminearum_oldlabel_negative']}`",
        "- positive definition: processed oldlabel positives reconstructed from the historical `gene_list.txt` replay audit and materialized under `data/processed/essential_gene/fgraminearum/oldlabel/positive_genes.tsv`",
        "- conceptual description: historical lethal + virulence old-label regime as preserved by the protocolized replay artifacts",
        f"- final counts: positives={old_row['positives']}, negatives={old_row['negatives']}, total={old_row['total']}",
        "",
        "## New Regime",
        f"- protocol: `{new_row['protocol_name']}`",
        "- role: intended mainline Fusarium benchmark regime",
        f"- positive source: `{config['label_sources']['fgraminearum_newlabel_positive']}`",
        f"- negative source: `{config['label_sources']['fgraminearum_newlabel_negative']}`",
        "- positive definition: processed newlabel positives = lethal PHI-supported positives union high-confidence yeast-transfer-supported positives",
        "- negative definition: processed newlabel negatives = weak-confidence none-derived genes after virulence and positive exclusions",
        "- conceptual description: current lethal + evolution newlabel regime materialized under `data/processed/essential_gene/fgraminearum/newlabel/`",
        f"- final counts: positives={new_row['positives']}, negatives={new_row['negatives']}, total={new_row['total']}",
        "",
        "## Unresolved Issues",
        "- Frozen protocol now consumes processed Fusarium label artifacts generated by `workflow/fgraminearum_label_materialization.smk`.",
        "- Reconstructing the newlabel negative-set provenance remains documented against the mirrored master evidence table under `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv`.",
    ]
    (labels_dir / "fgraminearum_label_protocol_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_split_protocol_summary(config: Dict[str, Any], split_rows: List[Dict[str, Any]]) -> None:
    splits_dir = Path(config["paths"]["splits_dir"])
    lines = [
        "# Frozen Split Protocol Summary",
        "",
        "- split type: stratified fixed split",
        f"- split seed: {int(config['runtime']['split_seed'])}",
        f"- test fraction: {float(config['runtime']['test_fraction']):.2f}",
        f"- val fraction: {float(config['runtime']['val_fraction']):.2f}",
        "- no model may generate its own split internally after this refactor",
        "",
        "| protocol | species | regime | train | val | test | train_pos | val_pos | test_pos | output |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in split_rows:
        lines.append(
            f"| {row['protocol_name']} | {row['species']} | {row['regime']} | {row['train']} | {row['val']} | {row['test']} | "
            f"{row['train_pos']} | {row['val_pos']} | {row['test_pos']} | `{row['output_file']}` |"
        )
    (splits_dir / "split_protocol_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_fusarium_split_summary(config: Dict[str, Any], split_rows: List[Dict[str, Any]]) -> None:
    splits_dir = Path(config["paths"]["splits_dir"])
    old_row = next(row for row in split_rows if row["protocol_name"] == "fgraminearum_oldlabel")
    new_row = next(row for row in split_rows if row["protocol_name"] == "fgraminearum_newlabel")
    lines = [
        "# Fusarium Split Protocol Summary",
        "",
        "## Shared Rules",
        "- split type: stratified fixed split",
        f"- split seed: {int(config['runtime']['split_seed'])}",
        f"- test fraction: {float(config['runtime']['test_fraction']):.2f}",
        f"- val fraction: {float(config['runtime']['val_fraction']):.2f}",
        "",
        "## Mainline",
        "- mainline split: `fgraminearum_newlabel_split.tsv`",
        f"- counts: train={new_row['train']}, val={new_row['val']}, test={new_row['test']}",
        f"- class balance: train_pos={new_row['train_pos']}, val_pos={new_row['val_pos']}, test_pos={new_row['test_pos']}",
        "",
        "## Legacy Comparison",
        "- legacy split: `fgraminearum_oldlabel_split.tsv`",
        f"- counts: train={old_row['train']}, val={old_row['val']}, test={old_row['test']}",
        f"- class balance: train_pos={old_row['train_pos']}, val_pos={old_row['val_pos']}, test_pos={old_row['test_pos']}",
        "",
        "## Caveats",
        "- The old regime exists only for replay and manuscript comparison.",
        "- The new regime is the only Fusarium split used by the new mainline benchmark protocol.",
    ]
    (splits_dir / "fgraminearum_split_protocol_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    labels_dir = Path(config["paths"]["labels_dir"])
    splits_dir = Path(config["paths"]["splits_dir"])
    labels_dir.mkdir(parents=True, exist_ok=True)
    splits_dir.mkdir(parents=True, exist_ok=True)

    label_rows: List[Dict[str, Any]] = []
    split_rows: List[Dict[str, Any]] = []

    for protocol_name in config["workflow"]["protocol_order"]:
        labels = build_label_manifest(config, protocol_name)
        label_path = label_output_path(config, protocol_name)
        labels.to_csv(label_path, sep="\t", index=False)
        label_counts = summarize_label_counts(labels)
        label_rows.append(
            {
                "protocol_name": protocol_name,
                "species": config["protocols"][protocol_name]["species"],
                "regime": config["protocols"][protocol_name]["regime"],
                "positives": label_counts["positives"],
                "negatives": label_counts["negatives"],
                "total": label_counts["total"],
                "source_manifest": labels["source_manifest"].iloc[0],
                "output_file": label_path.name,
                "is_mainline": bool(config["protocols"][protocol_name]["is_mainline"]),
            }
        )

        split_df = assign_splits(
            labels=labels,
            seed=int(config["runtime"]["split_seed"]),
            val_fraction=float(config["runtime"]["val_fraction"]),
            test_fraction=float(config["runtime"]["test_fraction"]),
        )
        split_path = split_output_path(config, protocol_name)
        split_df.to_csv(split_path, sep="\t", index=False)
        split_counts = summarize_split_counts(split_df)
        split_rows.append(
            {
                "protocol_name": protocol_name,
                "species": config["protocols"][protocol_name]["species"],
                "regime": config["protocols"][protocol_name]["regime"],
                "output_file": split_path.name,
                **split_counts,
            }
        )

    write_label_protocol_summary(config, label_rows)
    write_fusarium_label_summary(config, label_rows)
    write_split_protocol_summary(config, split_rows)
    write_fusarium_split_summary(config, split_rows)

    print(f"Wrote frozen labels to: {labels_dir}")
    print(f"Wrote frozen splits to: {splits_dir}")


if __name__ == "__main__":
    main()
