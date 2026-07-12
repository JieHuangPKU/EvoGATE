import argparse
from pathlib import Path

import pandas as pd

from src.eval.aggregate_frozen_protocol_runs import dataframe_to_markdown


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate per-run gate statistics from a figure output root")
    parser.add_argument("--output-root", required=True, type=str)
    parser.add_argument("--summary-dir", required=True, type=str)
    parser.add_argument("--prefix", required=True, type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    summary_dir = Path(args.summary_dir).resolve()
    summary_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for gate_path in sorted(output_root.glob("**/gate_statistics.tsv")):
        run_dir = gate_path.parent
        metrics_path = run_dir / "metrics.tsv"
        if not metrics_path.exists():
            continue
        metrics_row = pd.read_csv(metrics_path, sep="\t").iloc[0].to_dict()
        gate_df = pd.read_csv(gate_path, sep="\t")
        for row in gate_df.to_dict(orient="records"):
            rows.append(
                {
                    "protocol": metrics_row.get("protocol", ""),
                    "species": metrics_row.get("species", ""),
                    "regime": metrics_row.get("regime", ""),
                    "model": metrics_row.get("model", ""),
                    "feature_setting": metrics_row.get("feature_setting", ""),
                    "run_id": metrics_row.get("run_id", ""),
                    "seed": metrics_row.get("seed", ""),
                    **row,
                    "gate_statistics_path": str(gate_path),
                }
            )

    gate_stats = pd.DataFrame(rows)
    gate_tsv = summary_dir / f"{args.prefix}_gate_statistics.tsv"
    gate_md = summary_dir / f"{args.prefix}_gate_statistics.md"
    gate_stats.to_csv(gate_tsv, sep="\t", index=False)
    gate_md.write_text(
        dataframe_to_markdown(gate_stats, f"{args.prefix} Gate Statistics", "Per-run gate statistics exported by the gated GraphSAGE runs."),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
