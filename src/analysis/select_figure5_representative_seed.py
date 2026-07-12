import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


DEFAULT_PROTOCOLS = ["fgraminearum_newlabel", "scerevisiae"]
PROTOCOL_OUTPUT_SLUGS = {
    "fgraminearum_newlabel": "fgraminearum",
    "scerevisiae": "scerevisiae",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Select a representative Figure5 seed from five paired Bio vs Bio+ESM2 runs")
    parser.add_argument("--runtime-config", default="results/Figure3a/runtime/Figure3a_runtime_config.yaml", type=str)
    parser.add_argument("--upstream-root", default="outputs/Figure3a", type=str)
    parser.add_argument("--protocols", nargs="+", default=DEFAULT_PROTOCOLS)
    parser.add_argument("--subset", default="test", type=str)
    parser.add_argument("--output-table", required=True, type=str)
    parser.add_argument("--output-summary", required=True, type=str)
    return parser.parse_args()


def load_runtime_seed_list(runtime_config_path):
    config = yaml.safe_load(Path(runtime_config_path).read_text(encoding="utf-8"))
    return [int(seed) for seed in config["runtime"]["seed_list"]]


def protocol_output_slug(protocol):
    return PROTOCOL_OUTPUT_SLUGS.get(protocol, str(protocol))


def write_markdown(path, lines):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def mad_scale(values):
    values = np.asarray(values, dtype=float)
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    scale = mad if mad > 1e-12 else max(float(np.std(values, ddof=0)), 1.0)
    return median, scale


def load_metrics_row(upstream_root, protocol, feature_setting, seed):
    path = Path(upstream_root) / protocol / "GraphSAGE" / feature_setting / f"run_{seed}" / "metrics.tsv"
    return pd.read_csv(path, sep="\t").iloc[0].to_dict()


def load_prediction_df(upstream_root, protocol, feature_setting, seed, subset):
    path = Path(upstream_root) / protocol / "GraphSAGE" / feature_setting / f"run_{seed}" / "predictions.tsv"
    df = pd.read_csv(path, sep="\t")
    df = df[df["split"].astype(str) == subset].copy()
    df["graph_gene_id"] = df["graph_gene_id"].astype(str)
    df["label"] = df["label"].astype(float).astype(int)
    df["pred_label"] = df["pred_label"].astype(int)
    df["pred_score"] = df["pred_score"].astype(float)
    return df


def seed_metric_row(upstream_root, protocol, seed, subset):
    base_df = load_prediction_df(upstream_root, protocol, "ORT_EXP_SUB", seed, subset)
    esm2_df = load_prediction_df(upstream_root, protocol, "ORT_EXP_SUB_ESM2", seed, subset)
    merged = base_df.merge(
        esm2_df[["graph_gene_id", "pred_score", "pred_label"]],
        on="graph_gene_id",
        how="inner",
        validate="one_to_one",
        suffixes=("_baseline", "_esm2"),
    )
    rescued_mask = (
        (merged["label"].astype(int) == 1)
        & (merged["pred_label_baseline"].astype(int) == 0)
        & (merged["pred_label_esm2"].astype(int) == 1)
    )
    delta_scores = merged["pred_score_esm2"].astype(float) - merged["pred_score_baseline"].astype(float)
    esm2_metrics_row = load_metrics_row(upstream_root, protocol, "ORT_EXP_SUB_ESM2", seed)
    return {
        "rescued_gene_count": int(rescued_mask.sum()),
        "rescued_gene_fraction": float(rescued_mask.mean()),
        "rescued_mean_delta_probability": float(delta_scores.loc[rescued_mask].mean()) if rescued_mask.any() else 0.0,
        "all_gene_mean_delta_probability": float(delta_scores.mean()),
        "test_auprc_esm2": float(esm2_metrics_row["test_auprc"]),
        "test_mcc_esm2": float(esm2_metrics_row["test_mcc"]),
        "test_f1_esm2": float(esm2_metrics_row["test_f1"]),
        "test_auroc_esm2": float(esm2_metrics_row["test_auroc"]),
    }


def main():
    args = parse_args()
    seeds = load_runtime_seed_list(args.runtime_config)
    base_metric_names = [
        "rescued_gene_count",
        "rescued_gene_fraction",
        "rescued_mean_delta_probability",
        "all_gene_mean_delta_probability",
        "test_auprc_esm2",
        "test_mcc_esm2",
        "test_f1_esm2",
        "test_auroc_esm2",
    ]

    rows = []
    for seed in seeds:
        row = {"seed": seed}
        for protocol in args.protocols:
            prefix = protocol_output_slug(protocol)
            metrics_row = seed_metric_row(args.upstream_root, protocol, seed, args.subset)
            for metric_name, metric_value in metrics_row.items():
                row[f"{prefix}__{metric_name}"] = metric_value
        rows.append(row)

    selection_df = pd.DataFrame(rows).sort_values("seed").reset_index(drop=True)
    metric_columns = [column for column in selection_df.columns if column != "seed"]
    component_scores = []
    for metric in metric_columns:
        median, scale = mad_scale(selection_df[metric].to_numpy(dtype=float))
        deviation = (selection_df[metric].astype(float) - median).abs() / scale
        selection_df[f"{metric}_median"] = median
        selection_df[f"{metric}_scaled_abs_deviation"] = deviation
        component_scores.append(deviation)
    selection_df["aggregate_median_deviation_score"] = pd.concat(component_scores, axis=1).sum(axis=1)
    representative_idx = selection_df.sort_values(["aggregate_median_deviation_score", "seed"]).index[0]
    selection_df["is_representative_seed"] = False
    selection_df.loc[representative_idx, "is_representative_seed"] = True
    representative_seed = int(selection_df.loc[representative_idx, "seed"])

    output_table = Path(args.output_table)
    output_table.parent.mkdir(parents=True, exist_ok=True)
    selection_df.to_csv(output_table, sep="\t", index=False)

    representative_row = selection_df.loc[representative_idx]
    summary_lines = [
        "# Figure5 Representative Seed Selection",
        "",
        f"- Seeds considered: `{', '.join(str(seed) for seed in seeds)}`.",
        f"- Protocols considered: `{', '.join(args.protocols)}`.",
        f"- Base metric families: `{', '.join(base_metric_names)}`.",
        "- Species-specific metric columns are scored separately, then summed into a single cross-species deviation score.",
        f"- Selected representative seed: `{representative_seed}`.",
        "- Selection principle: closest overall behavior to the five-seed median, not the visually best-looking seed.",
        "- Exact rule: for every protocol-specific selection metric, compute the scaled absolute deviation from the five-seed median and choose the seed with the smallest aggregate cross-species deviation score.",
        f"- Representative seed score: `{representative_row['aggregate_median_deviation_score']:.6f}`.",
        f"- Full per-seed table: `{output_table}`.",
        "",
        "## Representative seed metrics",
        "",
    ]
    for protocol in args.protocols:
        prefix = protocol_output_slug(protocol)
        summary_lines.extend(
            [
                f"- `{prefix}` rescued gene count: `{int(representative_row[f'{prefix}__rescued_gene_count'])}`.",
                f"- `{prefix}` rescued gene fraction: `{representative_row[f'{prefix}__rescued_gene_fraction']:.6f}`.",
                f"- `{prefix}` mean rescued delta probability: `{representative_row[f'{prefix}__rescued_mean_delta_probability']:.6f}`.",
                f"- `{prefix}` Bio+ESM2 test AUPRC: `{representative_row[f'{prefix}__test_auprc_esm2']:.6f}`.",
                f"- `{prefix}` Bio+ESM2 test MCC: `{representative_row[f'{prefix}__test_mcc_esm2']:.6f}`.",
            ]
        )
    write_markdown(args.output_summary, summary_lines)


if __name__ == "__main__":
    main()
