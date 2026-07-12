import argparse
from pathlib import Path

import pandas as pd
import yaml

from src.eval.publication_summary import build_publication_summary


METHOD_ORDER = ["MLP", "RF", "SVM", "NB", "N2V_MLP", "DC", "CC"]
FEATURE_ORDER = ["ORT", "EXP", "SUB", "ORT_EXP", "ORT_SUB", "EXP_SUB", "ORT_EXP_SUB", "N2V", "network"]
METRIC_COLUMNS = ["auroc", "auprc", "mcc", "f1", "accuracy", "specificity"]


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate classical baseline benchmark outputs")
    parser.add_argument("--run-root", type=str)
    parser.add_argument("--aggregated-output", type=str)
    parser.add_argument("--output-root", type=str)
    parser.add_argument("--final-summary-output", type=str)
    return parser.parse_args()


def _load_yaml(path):
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _load_reference(run_dir):
    resolved_path = run_dir / "resolved_config.yaml"
    feature_summary_path = run_dir / "feature_summary.tsv"
    method_summary_path = run_dir / "method_summary.tsv"
    if resolved_path.exists():
        payload = _load_yaml(resolved_path)
        return {
            "species": payload["species"],
            "method": payload["method"],
            "feature_setting": payload["feature_setting"],
            "label_regime": payload["label_regime"],
        }
    if feature_summary_path.exists():
        row = pd.read_csv(feature_summary_path, sep="\t").iloc[0].to_dict()
        return {
            "species": row["species"],
            "method": row["method"],
            "feature_setting": row["feature_setting"],
            "label_regime": row["label_regime"],
        }
    if method_summary_path.exists():
        row = pd.read_csv(method_summary_path, sep="\t").iloc[0].to_dict()
        return {
            "species": row["species"],
            "method": row["method"],
            "feature_setting": row["feature_setting"],
            "label_regime": row["label_regime"],
        }
    raise FileNotFoundError(f"No reference metadata found under {run_dir}")


def _aggregate_one(run_root, aggregated_output):
    if not run_root.exists():
        raise FileNotFoundError(f"Run root does not exist: {run_root}")

    run_dirs = sorted(path for path in run_root.iterdir() if path.is_dir() and path.name.startswith("run_"))
    if run_dirs:
        metrics_frames = []
        reference = None
        for run_dir in run_dirs:
            metrics_path = run_dir / "metrics.tsv"
            if not metrics_path.exists():
                raise FileNotFoundError(f"Missing metrics.tsv under {run_dir}")
            metrics_frames.append(pd.read_csv(metrics_path, sep="\t"))
            if reference is None:
                reference = _load_reference(run_dir)
        metric_df = pd.concat(metrics_frames, ignore_index=True)
        out = dict(reference)
        out["n_runs"] = int(len(metric_df))
        for metric in METRIC_COLUMNS:
            out[f"{metric}_mean"] = float(metric_df[metric].mean())
            out[f"{metric}_std"] = float(metric_df[metric].std(ddof=0))
    else:
        metrics_path = run_root / "metrics.tsv"
        if not metrics_path.exists():
            raise FileNotFoundError(f"Missing metrics.tsv under {run_root}")
        metric_df = pd.read_csv(metrics_path, sep="\t")
        reference = _load_reference(run_root)
        row = metric_df.iloc[0].to_dict()
        out = dict(reference)
        out["n_runs"] = 1
        for metric in METRIC_COLUMNS:
            out[f"{metric}_mean"] = float(row[metric])
            out[f"{metric}_std"] = 0.0

    aggregated_output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([out]).to_csv(aggregated_output, sep="\t", index=False)


def _aggregate_final(output_root, final_summary_output):
    rows = []
    for species_dir in sorted(path for path in output_root.iterdir() if path.is_dir()):
        for method_dir in sorted(path for path in species_dir.iterdir() if path.is_dir()):
            for feature_dir in sorted(path for path in method_dir.iterdir() if path.is_dir()):
                agg_path = feature_dir / "aggregated_metrics.tsv"
                if not agg_path.exists():
                    continue
                rows.append(pd.read_csv(agg_path, sep="\t").iloc[0].to_dict())
    if rows:
        summary = pd.DataFrame(rows)
        summary["method_rank"] = summary["method"].map({name: idx for idx, name in enumerate(METHOD_ORDER)})
        summary["feature_rank"] = summary["feature_setting"].map({name: idx for idx, name in enumerate(FEATURE_ORDER)})
        summary = summary.sort_values(["species", "method_rank", "feature_rank"]).drop(columns=["method_rank", "feature_rank"]).reset_index(drop=True)
    else:
        summary = pd.DataFrame(
            columns=["species", "method", "feature_setting", "label_regime", "n_runs"] + [f"{metric}_{suffix}" for metric in METRIC_COLUMNS for suffix in ["mean", "std"]]
        )
    final_summary_output.parent.mkdir(parents=True, exist_ok=True)
    rename_map = {"method": "model"}
    for metric in METRIC_COLUMNS:
        rename_map[f"{metric}_mean"] = f"test_{metric}_mean"
        rename_map[f"{metric}_std"] = f"test_{metric}_std"
    publication_ready = summary.rename(columns=rename_map)
    build_publication_summary(publication_ready).to_csv(final_summary_output, sep="\t", index=False)


def main():
    args = parse_args()
    if args.run_root and args.aggregated_output:
        _aggregate_one(Path(args.run_root), Path(args.aggregated_output))
        return
    if args.output_root and args.final_summary_output:
        _aggregate_final(Path(args.output_root), Path(args.final_summary_output))
        return
    raise SystemExit("Provide either --run-root/--aggregated-output or --output-root/--final-summary-output")


if __name__ == "__main__":
    main()
