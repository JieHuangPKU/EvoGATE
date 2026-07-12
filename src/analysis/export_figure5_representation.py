import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from src.analysis.export_figure3d_representation import (
    DEFAULT_FEATURE_SETTINGS,
    DEFAULT_PROTOCOLS,
    DEFAULT_UPSTREAM_ROOTS,
    compute_embedding_metrics,
    fit_tsne_coords,
    fit_umap_coords,
    forward_penultimate_from_checkpoint,
    json_ready,
    load_protocol_dataset,
    load_run_context,
    locate_run_dir,
    normalize_model_name,
    output_stem,
    protocol_label_regime,
    resolve_graph_model_config,
    subset_indices,
    torch,
    yaml,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Export Figure5 representation-space plots and numeric tables")
    parser.add_argument("--base-config", default="configs/frozen_protocol.yaml", type=str)
    parser.add_argument("--protocols", nargs="+", default=DEFAULT_PROTOCOLS)
    parser.add_argument("--feature-settings", nargs="+", default=DEFAULT_FEATURE_SETTINGS)
    parser.add_argument("--model", default="GraphSAGE", type=str)
    parser.add_argument("--seed", default=1029, type=int)
    parser.add_argument("--subset", default="test", choices=["test", "val", "train", "all_labeled"])
    parser.add_argument("--plot-dir", default="results/Figure5/plots", type=str)
    parser.add_argument("--data-dir", default="results/Figure5/data", type=str)
    parser.add_argument("--table-dir", default="results/Figure5/tables", type=str)
    parser.add_argument("--summary-dir", default="results/Figure5/summary", type=str)
    parser.add_argument("--upstream-roots", nargs="+", default=DEFAULT_UPSTREAM_ROOTS)
    return parser.parse_args()


def save_plot(coords_df, x_col, y_col, projection_name, pdf_path, png_path):
    fig, ax = plt.subplots(figsize=(4.2, 4.0), facecolor="white")
    ax.set_facecolor("white")
    for label_value, label_name, color in [(0, "non-essential", "#4C78A8"), (1, "essential", "#E45756")]:
        subset = coords_df[coords_df["label"] == label_value].copy()
        ax.scatter(
            subset[x_col],
            subset[y_col],
            s=10,
            alpha=0.8,
            c=color,
            edgecolors="none",
            label=label_name,
        )
    ax.set_xlabel(f"{projection_name} 1")
    ax.set_ylabel(f"{projection_name} 2")
    ax.legend(frameon=False, loc="best")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path, format="pdf", dpi=300, transparent=False, bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=300, transparent=False, bbox_inches="tight")
    plt.close(fig)


def save_coords_table(subset_df, coords, output_path, protocol, feature_setting, model, seed, projection):
    out = pd.DataFrame(
        {
            "node_id": subset_df["graph_gene_id"].astype(str).to_numpy(),
            "species": subset_df["species"].astype(str).to_numpy(),
            "protocol": protocol,
            "label": subset_df["label"].astype(int).to_numpy(),
            "split": subset_df["split"].astype(str).to_numpy(),
            "feature_setting": feature_setting,
            "model": model,
            "seed": seed,
            "projection": projection,
            "dim1": coords[:, 0],
            "dim2": coords[:, 1],
        }
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False)
    return out


def main():
    args = parse_args()
    plot_dir = Path(args.plot_dir).resolve()
    data_dir = Path(args.data_dir).resolve()
    table_dir = Path(args.table_dir).resolve()
    summary_dir = Path(args.summary_dir).resolve()
    for path in [plot_dir, data_dir, table_dir, summary_dir]:
        path.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    metric_rows = []
    pairwise_inputs = {}

    for protocol in args.protocols:
        for feature_setting in args.feature_settings:
            run_dir = locate_run_dir(protocol, args.model, feature_setting, args.seed, args.upstream_roots)
            _, config_path = load_run_context(run_dir, args.base_config)
            base_config = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
            bundle = load_protocol_dataset(config_path, protocol, feature_setting)
            checkpoint = torch.load(run_dir / "best_model.pt", map_location="cpu")
            model_cfg = resolve_graph_model_config(base_config, normalize_model_name(args.model))
            x = torch.as_tensor(bundle["feature_matrix"], dtype=torch.float32)
            edge_index = torch.as_tensor(bundle["edge_index"].T, dtype=torch.long).contiguous()
            hidden, _ = forward_penultimate_from_checkpoint(
                x,
                edge_index,
                checkpoint["state_dict"],
                str(model_cfg.get("aggregator_type", "mean")).strip().lower(),
            )

            selected_idx = subset_indices(bundle, args.subset)
            subset_df = bundle["node_manifest"].iloc[selected_idx].copy().reset_index(drop=True)
            subset_df["species"] = bundle["species"]
            subset_df["label"] = subset_df["label"].astype(int)
            subset_df["split"] = subset_df["split"].astype(str)
            embedding = hidden.detach().cpu().numpy()[selected_idx]
            labels = subset_df["label"].to_numpy(dtype=int)

            projections = {
                "umap": fit_umap_coords(embedding),
                "tsne": fit_tsne_coords(embedding),
            }
            high_dim_metrics = compute_embedding_metrics(embedding, labels)
            stem = output_stem(protocol, feature_setting, args.subset, args.seed, args.model)

            pairwise_inputs.setdefault(protocol, {})[feature_setting] = {
                "split_manifest_path": bundle["split_manifest_path"],
                "label_manifest_path": bundle["label_manifest_path"],
                "subset": args.subset,
                "node_ids": subset_df["graph_gene_id"].astype(str).tolist(),
                "labels": subset_df["label"].astype(int).tolist(),
            }

            row = {
                "protocol": protocol,
                "species": bundle["species"],
                "label_regime": protocol_label_regime(protocol),
                "feature_setting": feature_setting,
                "model": args.model,
                "seed": int(args.seed),
                "subset": args.subset,
                "node_count": int(len(subset_df)),
                "checkpoint_path": str((run_dir / "best_model.pt").resolve()),
                "split_manifest_path": bundle["split_manifest_path"],
                "label_manifest_path": bundle["label_manifest_path"],
            }
            for projection, (coords, params) in projections.items():
                coords_path = data_dir / f"Figure5_representation_{stem}_{projection}_coords.tsv"
                pdf_path = plot_dir / f"Figure5_representation_{stem}_{projection}.pdf"
                png_path = plot_dir / f"Figure5_representation_{stem}_{projection}.png"
                coords_df = save_coords_table(subset_df, coords, coords_path, protocol, feature_setting, args.model, args.seed, projection)
                save_plot(coords_df, "dim1", "dim2", projection.upper() if projection == "umap" else "t-SNE", pdf_path, png_path)
                metrics = compute_embedding_metrics(coords, labels)
                metric_rows.extend(
                    [
                        {
                            "protocol": protocol,
                            "species": bundle["species"],
                            "feature_setting": feature_setting,
                            "model": args.model,
                            "seed": args.seed,
                            "subset": args.subset,
                            "projection": projection,
                            "metric": metric_name,
                            "value": metric_value,
                            "n": 1,
                            "mean": metric_value,
                            "std": 0.0,
                            "variance": 0.0,
                            "coords_tsv": str(coords_path),
                            "plot_pdf": str(pdf_path),
                            "plot_png": str(png_path),
                        }
                        for metric_name, metric_value in metrics.items()
                    ]
                )
                row[f"{projection}_coords_tsv"] = str(coords_path)
                row[f"{projection}_plot_pdf"] = str(pdf_path)
                row[f"{projection}_plot_png"] = str(png_path)
                row[f"{projection}_params_json"] = json.dumps(json_ready(params), sort_keys=True)

            for metric_name, metric_value in high_dim_metrics.items():
                metric_rows.append(
                    {
                        "protocol": protocol,
                        "species": bundle["species"],
                        "feature_setting": feature_setting,
                        "model": args.model,
                        "seed": args.seed,
                        "subset": args.subset,
                        "projection": "hidden_high_dim",
                        "metric": metric_name,
                        "value": metric_value,
                        "n": 1,
                        "mean": metric_value,
                        "std": 0.0,
                        "variance": 0.0,
                        "coords_tsv": "",
                        "plot_pdf": "",
                        "plot_png": "",
                    }
                )
                row[f"hidden_high_dim_{metric_name}"] = metric_value

            manifest_rows.append(row)

    manifest_path = table_dir / "Figure5_representation_manifest.tsv"
    metrics_path = table_dir / "Figure5_representation_projection_metrics.tsv"
    pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)
    pd.DataFrame(metric_rows).to_csv(metrics_path, sep="\t", index=False)

    pairwise_lines = []
    for protocol, feature_map in sorted(pairwise_inputs.items()):
        if len(feature_map) < 2:
            continue
        baseline = feature_map[args.feature_settings[0]]
        compare = feature_map[args.feature_settings[1]]
        pairwise_lines.append(
            "- `{0}`: same split manifest=`{1}`, same label manifest=`{2}`, same node ids=`{3}`, same labels=`{4}`".format(
                protocol,
                baseline["split_manifest_path"] == compare["split_manifest_path"],
                baseline["label_manifest_path"] == compare["label_manifest_path"],
                baseline["node_ids"] == compare["node_ids"],
                baseline["labels"] == compare["labels"],
            )
        )

    summary_lines = [
        "# Figure5 Representation Export",
        "",
        "- Scope: single-species representation-space exports for the Figure5 mechanism module.",
        "- Model: `GraphSAGE`.",
        "- Feature settings: `ORT_EXP_SUB`, `ORT_EXP_SUB_ESM2`.",
        "- Output manifest: `{0}`".format(manifest_path),
        "- Projection metrics table: `{0}`".format(metrics_path),
        "",
        "## Pairwise alignment checks",
        "",
    ]
    summary_lines.extend(pairwise_lines or ["- No pairwise validation rows were generated."])
    (summary_dir / "Figure5_representation_export.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
