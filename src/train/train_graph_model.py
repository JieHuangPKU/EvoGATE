from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from src.features.load_embeddings import load_embedding_index, load_feature_matrix
from src.models.graph_models import build_graph_model, inspect_graph_dependencies
from src.schemas.graph_schema import validate_graph_manifest_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="First real graph prototype trainer for ProGATE_v2")
    parser.add_argument("--config", type=str, required=True, help="Path to graph prototype YAML config")
    return parser.parse_args()


def load_config(config_path: str | Path) -> dict[str, Any]:
    with Path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required input table not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _bool_mask(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def _prepare_edge_index(edges_df: pd.DataFrame, nodes_df: pd.DataFrame):
    import torch

    index_lookup = {canonical_id: idx for idx, canonical_id in enumerate(nodes_df["canonical_gene_id"].astype(str))}
    kept_edges = edges_df[
        edges_df["source_canonical_gene_id"].astype(str).isin(index_lookup)
        & edges_df["target_canonical_gene_id"].astype(str).isin(index_lookup)
    ].copy()

    src = kept_edges["source_canonical_gene_id"].astype(str).map(index_lookup).astype(int).to_numpy()
    dst = kept_edges["target_canonical_gene_id"].astype(str).map(index_lookup).astype(int).to_numpy()
    src_bidirectional = list(src) + list(dst)
    dst_bidirectional = list(dst) + list(src)
    edge_index = torch.tensor([src_bidirectional, dst_bidirectional], dtype=torch.long)
    return edge_index, kept_edges


def _write_inputs_summary(output_path: Path, rows: list[dict[str, Any]]) -> None:
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


def _seed_everything(seed: int) -> None:
    import random
    import numpy as np
    import torch

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    output_dir = Path(config["paths"]["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    graph_manifest_path = Path(config["paths"]["graph_manifest_path"])
    manifest_df = validate_graph_manifest_file(graph_manifest_path)
    graph_id = str(config["runtime"]["graph_id"])
    model_name = str(config["runtime"]["model_name"]).strip().lower()
    graph_rows = manifest_df[manifest_df["graph_id"].astype(str) == graph_id].copy()
    if graph_rows.empty:
        raise ValueError(f"graph_id not found in graph manifest: {graph_id}")
    graph_row = graph_rows.iloc[0]

    nodes_path = Path(graph_row["adapted_nodes_path"] or config["paths"]["graph_nodes_path"])
    edges_path = Path(graph_row["adapted_edges_path"] or config["paths"]["graph_edges_path"])
    unmatched_path = Path(graph_row["adapted_unmatched_nodes_path"] or config["paths"]["graph_unmatched_nodes_path"])

    print(f"[train_graph_model] graph load graph_id={graph_id}")
    nodes_df = _read_tsv(nodes_path)
    edges_df = _read_tsv(edges_path)
    unmatched_df = _read_tsv(unmatched_path)
    print(
        f"[train_graph_model] graph loaded nodes={len(nodes_df)} edges={len(edges_df)} unmatched={len(unmatched_df)}"
    )

    print("[train_graph_model] feature load start")
    embedding_index = load_embedding_index(
        config["paths"]["feature_manifest_path"],
        source_name=config["runtime"]["feature_source_name"],
    )
    feature_input = nodes_df[["species", "canonical_gene_id"]].copy()
    x_np, aligned_nodes = load_feature_matrix(
        feature_input,
        embedding_index,
        require_all=True,
        require_pooled_features=bool(config["runtime"].get("require_pooled_features", True)),
    )
    if list(aligned_nodes["canonical_gene_id"].astype(str)) != list(nodes_df["canonical_gene_id"].astype(str)):
        raise ValueError("Feature-aligned node order does not match adapted graph node order")
    print(f"[train_graph_model] feature load finish feature_dim={x_np.shape[1]}")

    import torch

    device_name = str(config["runtime"].get("device", "cpu")).strip().lower()
    device = torch.device("cuda" if device_name == "cuda" and torch.cuda.is_available() else "cpu")
    _seed_everything(int(config["runtime"]["seed"]))

    x = torch.tensor(x_np, dtype=torch.float32, device=device)
    edge_index, kept_edges = _prepare_edge_index(edges_df, nodes_df)
    edge_index = edge_index.to(device)

    benchmark_membership = nodes_df["benchmark_membership"].astype(str)
    has_feature_mask = torch.tensor(_bool_mask(nodes_df["has_feature"]).to_numpy(), dtype=torch.bool, device=device)
    has_label_mask = torch.tensor(_bool_mask(nodes_df["has_label"]).to_numpy(), dtype=torch.bool, device=device)
    broad_mask = torch.tensor(benchmark_membership.str.contains("broad79", regex=False).to_numpy(), dtype=torch.bool, device=device)
    strict_mask = torch.tensor(benchmark_membership.str.contains("strict29", regex=False).to_numpy(), dtype=torch.bool, device=device)
    conflict_mask = torch.tensor(benchmark_membership.str.contains("conflict8", regex=False).to_numpy(), dtype=torch.bool, device=device)
    inference_mask = has_feature_mask.clone()
    benchmark_only_mask = broad_mask | strict_mask | conflict_mask
    train_mask = torch.zeros(len(nodes_df), dtype=torch.bool, device=device)
    val_mask = torch.zeros(len(nodes_df), dtype=torch.bool, device=device)
    test_mask = torch.zeros(len(nodes_df), dtype=torch.bool, device=device)
    anchor_mask = has_label_mask.clone()

    print(
        "[train_graph_model] mask summary "
        f"inference={int(inference_mask.sum().item())} anchors={int(anchor_mask.sum().item())} "
        f"broad79={int(broad_mask.sum().item())} strict29={int(strict_mask.sum().item())} conflict8={int(conflict_mask.sum().item())}"
    )

    dependency_status = inspect_graph_dependencies()
    backend = str(config["runtime"].get("backend", "auto")).strip().lower()
    train_mode = "unsupervised_edge_reconstruction_plus_positive_anchor_scoring"
    model = build_graph_model(
        model_name,
        input_dim=int(x.shape[1]),
        hidden_dim=int(config["runtime"]["hidden_dim"]),
        output_dim=int(config["runtime"]["embedding_dim"]),
        num_layers=int(config["runtime"]["num_layers"]),
        dropout=float(config["runtime"]["dropout"]),
        backend=backend,
    ).to(device)
    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=float(config["runtime"]["lr"]),
        weight_decay=float(config["runtime"].get("weight_decay", 1e-5)),
    )

    epoch_rows: list[dict[str, Any]] = []
    epochs = int(config["runtime"]["epochs"])
    print(f"[train_graph_model] train start mode={train_mode} epochs={epochs} backend={backend}")
    start_time = time.time()
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()
        embeddings = model.encode(x, edge_index)
        loss = model.reconstruction_loss(embeddings, edge_index)
        loss.backward()
        optimizer.step()

        elapsed = time.time() - start_time
        epoch_rows.append(
            {
                "epoch": epoch,
                "loss": float(loss.detach().cpu().item()),
                "elapsed_sec": round(elapsed, 4),
            }
        )
        print(f"[train_graph_model] epoch={epoch}/{epochs} loss={loss.detach().cpu().item():.6f}")

    model.eval()
    with torch.no_grad():
        final_embeddings = model.encode(x, edge_index)
        pred_score = model.decode_anchor_scores(
            final_embeddings,
            anchor_mask=anchor_mask,
            temperature=float(config["runtime"].get("anchor_temperature", 5.0)),
        )
    pred_score_cpu = pred_score.detach().cpu().numpy()
    embedding_cpu = final_embeddings.detach().cpu().numpy()

    prediction_df = nodes_df.copy()
    prediction_df["graph_id"] = graph_id
    prediction_df["model_name"] = model_name
    prediction_df["backend"] = model.backend_name
    prediction_df["prediction_stage"] = train_mode
    prediction_df["pred_score"] = pred_score_cpu
    prediction_df = prediction_df.sort_values("pred_score", ascending=False, kind="stable").reset_index(drop=True)
    prediction_df["pred_rank"] = prediction_df.index + 1
    prediction_df["runtime_note"] = (
        "Scores come from unsupervised GraphSAGE embeddings ranked by positive-anchor affinity; "
        "no fake negatives were introduced."
    )

    model_path = output_dir / "graphsage_model.pt"
    torch.save(
        {
            "graph_id": graph_id,
            "model_name": model_name,
            "backend": model.backend_name,
            "state_dict": model.state_dict(),
            "model_summary": model.model_summary(),
            "config": config,
        },
        model_path,
    )
    pd.DataFrame(epoch_rows).to_csv(output_dir / "graphsage_runtime_metrics.tsv", sep="\t", index=False)
    prediction_df.to_csv(output_dir / "graphsage_predictions.tsv", sep="\t", index=False)
    prediction_df.to_csv(output_dir / "graph_predictions.tsv", sep="\t", index=False)
    (output_dir / "graphsage_node_embeddings.npy").write_bytes(b"")  # placeholder file for path stability
    # Save embeddings using numpy after the path exists.
    import numpy as np

    np.save(output_dir / "graphsage_node_embeddings.npy", embedding_cpu)

    inputs_rows = [
        {"item": "graph_id", "value": graph_id},
        {"item": "model_name", "value": model_name},
        {"item": "backend", "value": model.backend_name},
        {"item": "node_count", "value": len(nodes_df)},
        {"item": "edge_count", "value": len(kept_edges)},
        {"item": "bidirectional_edge_count", "value": int(edge_index.shape[1])},
        {"item": "feature_dim", "value": int(x.shape[1])},
        {"item": "feature_ready_nodes", "value": int(has_feature_mask.sum().item())},
        {"item": "gold_label_nodes", "value": int(has_label_mask.sum().item())},
        {"item": "anchor_nodes", "value": int(anchor_mask.sum().item())},
        {"item": "broad79_nodes", "value": int(broad_mask.sum().item())},
        {"item": "strict29_nodes", "value": int(strict_mask.sum().item())},
        {"item": "conflict8_nodes", "value": int(conflict_mask.sum().item())},
        {"item": "inference_mask_nodes", "value": int(inference_mask.sum().item())},
        {"item": "benchmark_only_mask_nodes", "value": int(benchmark_only_mask.sum().item())},
        {"item": "train_mask_nodes", "value": int(train_mask.sum().item())},
        {"item": "val_mask_nodes", "value": int(val_mask.sum().item())},
        {"item": "test_mask_nodes", "value": int(test_mask.sum().item())},
        {"item": "unmatched_baseline_nodes", "value": len(unmatched_df)},
        {"item": "train_mode", "value": train_mode},
        {"item": "epoch_count", "value": epochs},
        {"item": "final_loss", "value": epoch_rows[-1]["loss"]},
    ]
    _write_inputs_summary(output_dir / "graph_train_inputs_summary.tsv", inputs_rows)

    summary_lines = [
        "# GraphSAGE Train Summary",
        "",
        "## Run",
        f"- graph_id: {graph_id}",
        f"- model_name: {model_name}",
        f"- backend: {model.backend_name}",
        f"- device: {device.type}",
        f"- train_mode: {train_mode}",
        "",
        "## Loading",
        f"- nodes: {len(nodes_df)}",
        f"- edges: {len(kept_edges)}",
        f"- feature_dim: {int(x.shape[1])}",
        f"- unmatched_baseline_nodes: {len(unmatched_df)}",
        "",
        "## Masks",
        f"- inference_mask: {int(inference_mask.sum().item())}",
        f"- anchor_mask: {int(anchor_mask.sum().item())}",
        f"- broad79_mask: {int(broad_mask.sum().item())}",
        f"- strict29_mask: {int(strict_mask.sum().item())}",
        f"- conflict8_mask: {int(conflict_mask.sum().item())}",
        f"- train_mask: {int(train_mask.sum().item())}",
        f"- val_mask: {int(val_mask.sum().item())}",
        f"- test_mask: {int(test_mask.sum().item())}",
        "",
        "## Training",
        "- This prototype performs real parameter updates, but the updates are driven by unsupervised edge reconstruction rather than supervised essential-vs-nonessential labels.",
        f"- epochs: {epochs}",
        f"- final_loss: {epoch_rows[-1]['loss']}",
        f"- model_path: {model_path}",
        f"- predictions_path: {output_dir / 'graphsage_predictions.tsv'}",
        "",
        "## Interpretation",
        "- This is a real training prototype, not just forward-only scoring, but it is still inference-first and ranking-only on Fusarium.",
        "- No unresolved genes were converted into negatives.",
    ]
    (output_dir / "graphsage_train_summary.md").write_text("\n".join(summary_lines), encoding="utf-8")

    runtime_lines = [
        "# Graph Runtime Check",
        "",
        "## Run",
        f"- graph_id: {graph_id}",
        f"- model_name: {model_name}",
        f"- backend: {model.backend_name}",
        f"- stage: first real GraphSAGE prototype",
        f"- graph_heavy_ready: no",
        "",
        "## Dependency Check",
    ]
    for note in dependency_status.notes:
        runtime_lines.append(f"- {note}")
    runtime_lines.extend(
        [
            "",
            "## Outcome",
            f"- model_build: ok",
            f"- epoch_count: {epochs}",
            f"- final_loss: {epoch_rows[-1]['loss']}",
            f"- prediction_table: {output_dir / 'graphsage_predictions.tsv'}",
            "",
            "## Interpretation",
            "- The prototype completed real parameter updates and produced scored node predictions.",
            "- The training signal came from graph structure reconstruction, while the final node ranking came from positive-anchor affinity in the learned embedding space.",
        ]
    )
    (output_dir / "runtime_check.md").write_text("\n".join(runtime_lines), encoding="utf-8")

    print(f"[train_graph_model] model save {model_path}")
    print(f"[train_graph_model] prediction save {output_dir / 'graphsage_predictions.tsv'}")
    print(f"[train_graph_model] summary save {output_dir / 'graphsage_train_summary.md'}")


if __name__ == "__main__":
    main()
