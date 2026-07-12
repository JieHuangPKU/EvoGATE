from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from src.schemas.graph_schema import GENE_LEVEL_GRAPH_TYPES, validate_graph_manifest_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize graph-ready dataset compatibility for ProGATE_v2")
    parser.add_argument("--config", type=str, required=True, help="Path to graph-ready YAML config")
    return parser.parse_args()


def load_config(config_path: str | Path) -> dict[str, Any]:
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required input table not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    output_dir = Path(config["paths"]["graph_ready_output_dir"])
    manifest_df = validate_graph_manifest_file(output_dir / "graph_manifest.tsv")

    support_df = _read_tsv(Path(config["paths"]["baseline_dataset_dir"]) / "support_supervised_samples.tsv")
    inference_df = _read_tsv(Path(config["paths"]["baseline_dataset_dir"]) / "fgraminearum_inference_pool.tsv")

    combined_nodes = pd.concat([support_df, inference_df], ignore_index=True)
    species_node_counts = combined_nodes.groupby("species", sort=True)["canonical_gene_id"].nunique().to_dict()

    manifest_df["runtime_ready_bool"] = _bool_series(manifest_df["runtime_ready"])
    manifest_df["compatible_bool"] = _bool_series(manifest_df["compatible_with_baseline"])
    manifest_df["asset_exists_bool"] = _bool_series(manifest_df["asset_exists"])
    if "adapter_ready" in manifest_df.columns:
        manifest_df["adapter_ready_bool"] = _bool_series(manifest_df["adapter_ready"])
    else:
        manifest_df["adapter_ready_bool"] = False

    gene_ready_species: list[str] = []
    current_blockers: list[str] = []
    schema_only_rows: list[str] = []
    adapter_ready_rows: list[str] = []

    for species, species_df in manifest_df.groupby("species", sort=True):
        gene_placeholder = species_df[species_df["graph_type"] == "gene_graph"].copy()
        if not gene_placeholder.empty and gene_placeholder["baseline_node_count"].astype(int).iloc[0] > 0:
            gene_ready_species.append(
                f"{species}: node_universe={int(gene_placeholder['baseline_node_count'].iloc[0])}, "
                f"feature_ready={int(gene_placeholder['feature_ready_node_count'].iloc[0])}"
            )

        for _, row in species_df.iterrows():
            if row["graph_type"] == "gene_graph":
                continue
            if row["adapter_ready_bool"]:
                adapter_ready_rows.append(
                    f"{row['species']}::{row['graph_type']}: "
                    f"nodes={row.get('adapter_node_count', '0')}, edges={row.get('adapter_edge_count', '0')}"
                )
            label = f"{row['species']}::{row['graph_type']}"
            if not row["runtime_ready_bool"]:
                reason = row["missing_reason"] or "not runtime-ready yet"
                schema_only_rows.append(f"{label}: {reason}")

        ppi_rows = species_df[species_df["graph_type"] == "ppi_graph"].copy()
        if not ppi_rows.empty:
            row = ppi_rows.iloc[0]
            if row["asset_exists_bool"] and not row["compatible_bool"]:
                current_blockers.append(
                    f"{species} ppi_graph: asset exists but node ids are not yet canonical-id aligned"
                )

        orthology_rows = species_df[species_df["graph_type"] == "orthology_graph"].copy()
        if not orthology_rows.empty:
            row = orthology_rows.iloc[0]
            if row["asset_exists_bool"] and not row["compatible_bool"]:
                current_blockers.append(
                    f"{species} orthology_graph: asset exists but current matrix ids are not directly compatible with canonical_gene_id"
                )

    current_blockers.extend(
        [
            "residue_graph assets exist in Bingo raw `.pt` form, but there is no v2 residue-graph builder or residue-level runtime yet",
            "Fusarium still has 15 unresolved canonical gaps, so any future gene-level graph adapter must preserve the existing non-negative handling",
            "no graph-heavy trainer is implemented in v2 yet by design",
        ]
    )

    graph_status_df = manifest_df[
        [
            "species",
            "graph_type",
            "asset_exists",
            "compatible_with_baseline",
            "runtime_ready",
            "adapter_ready",
            "baseline_node_count",
            "feature_ready_node_count",
            "asset_overlap_count",
            "missing_reason",
        ]
    ].copy()
    graph_status_df.to_csv(output_dir / "graph_dataset_status.tsv", sep="\t", index=False)

    summary_lines = [
        "# Graph-Ready Dataset Summary",
        "",
        "## Direct Answers",
        f"1. graph-ready dataset summary generated: yes",
        f"2. species with gene-level graph-ready node universes: {', '.join(sorted(species_node_counts))}",
        "3. graph-heavy training entered: no",
        "",
        "## Gene-Level Graph-Ready Coverage",
    ]
    for line in gene_ready_species:
        summary_lines.append(f"- {line}")

    summary_lines.extend(["", "## Runtime-Ready Graph Rows"])
    ready_rows = manifest_df[manifest_df["runtime_ready_bool"]].copy()
    if ready_rows.empty:
        summary_lines.append("- none")
    else:
        for _, row in ready_rows.iterrows():
            summary_lines.append(
                f"- {row['species']}::{row['graph_type']}: "
                f"asset_exists={row['asset_exists']}, compatible_with_baseline={row['compatible_with_baseline']}"
            )

    summary_lines.extend(["", "## Adapter-Ready Gene Graph Rows"])
    if not adapter_ready_rows:
        summary_lines.append("- none")
    else:
        for line in adapter_ready_rows:
            summary_lines.append(f"- {line}")

    summary_lines.extend(["", "## Schema-Ready But Not Runtime-Ready"])
    for line in schema_only_rows:
        summary_lines.append(f"- {line}")

    summary_lines.extend(["", "## Current Blockers"])
    for line in current_blockers:
        summary_lines.append(f"- {line}")

    summary_lines.extend(
        [
            "",
            "## Interpretation",
            "- All four species already have canonical-gene node universes derived from the baseline dataset, so a gene-level graph-ready dataset can be described now without changing the current dataset contract.",
            "- Adapter-ready means a raw gene-level graph asset has already been canonicalized into explicit node/edge tables under `outputs/graph_ready/`.",
            "- `ppi_graph` and `orthology_graph` are the closest future graph types, but only the Fusarium candidate assets are already close to canonical-id compatibility; the others still need a mapping layer.",
            "- `residue_graph` is intentionally separate: the raw Bingo `.pt` assets are useful future inputs, but they must not be mixed into the gene-level graph runtime defined in this round.",
        ]
    )

    (output_dir / "graph_dataset_summary.md").write_text("\n".join(summary_lines), encoding="utf-8")
    (output_dir / "graph_dataset_summary.json").write_text(
        json.dumps(
            {
                "gene_ready_species": gene_ready_species,
                "schema_only_rows": schema_only_rows,
                "current_blockers": current_blockers,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    print(f"Wrote graph dataset summary to: {output_dir / 'graph_dataset_summary.md'}")
    print(f"Wrote graph dataset status table to: {output_dir / 'graph_dataset_status.tsv'}")


if __name__ == "__main__":
    main()
