import argparse
from pathlib import Path

from src.eval.aggregate_frozen_protocol_runs import aggregate_runs, collect_metrics, dataframe_to_markdown
from src.eval.publication_summary import build_publication_summary, publication_markdown


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Figure-prefixed benchmark summary tables from a frozen-protocol output root")
    parser.add_argument("--output-root", required=True, type=str)
    parser.add_argument("--summary-dir", required=True, type=str)
    parser.add_argument("--prefix", required=True, type=str)
    parser.add_argument("--target-name", required=True, type=str)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root).resolve()
    summary_dir = Path(args.summary_dir).resolve()
    summary_dir.mkdir(parents=True, exist_ok=True)

    per_run = collect_metrics(output_root)
    if "target" not in per_run.columns:
        per_run["target"] = args.target_name
    else:
        per_run["target"] = per_run["target"].fillna(args.target_name)

    aggregated = aggregate_runs(per_run)
    if "target" not in aggregated.columns:
        aggregated["target"] = args.target_name
    else:
        aggregated["target"] = aggregated["target"].fillna(args.target_name)

    publication = build_publication_summary(aggregated)

    per_run_tsv = summary_dir / f"{args.prefix}_per_run_metrics.tsv"
    aggregated_tsv = summary_dir / f"{args.prefix}_aggregated_metrics.tsv"
    final_tsv = summary_dir / f"{args.prefix}_final_summary.tsv"
    per_run_md = summary_dir / f"{args.prefix}_per_run_metrics.md"
    aggregated_md = summary_dir / f"{args.prefix}_aggregated_metrics.md"
    final_md = summary_dir / f"{args.prefix}_final_summary.md"

    per_run.to_csv(per_run_tsv, sep="\t", index=False)
    aggregated.to_csv(aggregated_tsv, sep="\t", index=False)
    publication.to_csv(final_tsv, sep="\t", index=False)

    per_run_md.write_text(
        dataframe_to_markdown(per_run, f"{args.prefix} Per-Run Metrics", "Per-run benchmark rows collected from the figure-specific output root."),
        encoding="utf-8",
    )
    aggregated_md.write_text(
        dataframe_to_markdown(aggregated, f"{args.prefix} Aggregated Metrics", "Aggregated benchmark rows grouped by the frozen-protocol benchmark identity columns."),
        encoding="utf-8",
    )
    final_md.write_text(
        publication_markdown(
            publication,
            f"{args.prefix} Final Summary",
            "Publication-facing benchmark summary with fixed split, fixed seeds, and mean ± std test metrics.",
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
