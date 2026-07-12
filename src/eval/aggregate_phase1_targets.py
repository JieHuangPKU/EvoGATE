from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.eval.publication_summary import build_publication_summary, publication_markdown


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate per-target Phase 1 benchmark summaries")
    parser.add_argument("--results-root", required=True, help="Legacy phase1 benchmark results root")
    parser.add_argument("--output-dir", required=True, help="Root summary directory")
    return parser.parse_args()


def _collect_tables(results_root: Path, filename: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for target_dir in sorted(path for path in results_root.iterdir() if path.is_dir()):
        table_path = target_dir / filename
        if table_path.exists():
            frames.append(pd.read_csv(table_path, sep="\t"))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    args = parse_args()
    results_root = Path(args.results_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    per_run = _collect_tables(results_root, "per_run_metrics.tsv")
    aggregated = _collect_tables(results_root, "aggregated_metrics.tsv")
    publication = build_publication_summary(aggregated)

    per_run.to_csv(output_dir / "per_run_metrics.tsv", sep="\t", index=False)
    aggregated.to_csv(output_dir / "aggregated_metrics.tsv", sep="\t", index=False)
    publication.to_csv(output_dir / "final_summary.tsv", sep="\t", index=False)
    (output_dir / "per_run_metrics.md").write_text(
        "# Per-Run Metrics\n\n" + (per_run.to_markdown(index=False) if not per_run.empty else "No per-run rows available.") + "\n",
        encoding="utf-8",
    )
    (output_dir / "aggregated_metrics.md").write_text(
        "# Aggregated Metrics\n\n" + (aggregated.to_markdown(index=False) if not aggregated.empty else "No aggregated rows available.") + "\n",
        encoding="utf-8",
    )
    (output_dir / "final_summary.md").write_text(
        publication_markdown(
            publication,
            "Final Summary",
            "Publication-facing Phase 1 benchmark summary across all six benchmark targets under the frozen split contract.",
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
