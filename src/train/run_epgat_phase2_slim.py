"""
Run the slim Phase 2 matrix:
- species: human, celegans, fgraminearum_canonical
- feature_mode: baseline, baseline_plus_esm2, baseline_plus_prott5
"""

import argparse
import os
import subprocess

import pandas as pd
import yaml


CONFIGS = [
    "configs/epgat_extended_human.yaml",
    "configs/epgat_extended_celegans.yaml",
    "configs/epgat_extended_fgraminearum_canonical.yaml",
]
MODES = ["baseline", "baseline_plus_esm2", "baseline_plus_prott5"]


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def main():
    parser = argparse.ArgumentParser(description="Run slim Phase 2 extended matrix")
    parser.add_argument("--output-summary", default="outputs/epgat_extended/phase2_slim_summary.tsv", type=str)
    args = parser.parse_args()

    os.makedirs("outputs/epgat_extended", exist_ok=True)
    summary_rows = []

    baseline_lookup = {}

    for config_path in CONFIGS:
        config = load_yaml(config_path)
        species = config["species"]["name"]
        for mode in MODES:
            subprocess.check_call(["python", "src/data/build_epgat_extended_dataset.py", "--config", config_path, "--feature-mode", mode])
            subprocess.check_call(["python", "src/train/train_epgat_extended.py", "--config", config_path, "--feature-mode", mode])
            subprocess.check_call(["python", "src/eval/evaluate_epgat_extended.py", "--config", config_path, "--feature-mode", mode])

            run_dir = os.path.join(config["paths"]["output_root"], "{}__{}".format(config["run"]["name"], mode))
            metrics = pd.read_csv(os.path.join(run_dir, "metrics.tsv"), sep="\t").iloc[0]
            coverage = pd.read_csv(os.path.join(run_dir, "plm_coverage.tsv"), sep="\t").iloc[0]
            audit = pd.read_csv(os.path.join(run_dir, "dataset_alignment_audit.tsv"), sep="\t").iloc[0]

            row = {
                "species": species,
                "feature_mode": mode,
                "num_nodes_final": int(audit["row_count_after_join"]) if "row_count_after_join" in audit else int(audit["canonical_nodes_final"]),
                "num_labeled": int((pd.read_csv(os.path.join(run_dir, "label_manifest.tsv"), sep="\t")["is_labeled"].astype(str).str.lower().isin(["true","1","yes"])).sum()),
                "plm_name": coverage["plm_name"],
                "plm_dim": int(coverage["plm_dim"]),
                "plm_coverage": float(coverage["coverage"]),
                "auroc": float(metrics["auroc"]),
                "auprc": float(metrics["auprc"]),
                "accuracy": float(metrics["accuracy"]),
                "f1": float(metrics["f1"]),
                "mcc": float(metrics["mcc"]),
                "delta_vs_baseline_auroc": 0.0,
                "delta_vs_baseline_auprc": 0.0,
                "delta_vs_baseline_mcc": 0.0,
                "notes": "",
            }
            if mode == "baseline":
                baseline_lookup[species] = row
            else:
                base = baseline_lookup[species]
                row["delta_vs_baseline_auroc"] = row["auroc"] - base["auroc"]
                row["delta_vs_baseline_auprc"] = row["auprc"] - base["auprc"]
                row["delta_vs_baseline_mcc"] = row["mcc"] - base["mcc"]
            summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(args.output_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
