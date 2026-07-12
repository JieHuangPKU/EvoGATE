# Classical Baseline Migration

## Scope

This migration adds the non-GNN classical EPGAT baselines to `ProGATE_v2` as a reproducible Snakemake workflow.

Migrated methods:

- `MLP`
- `RF`
- `SVM`
- `NB`
- `N2V_MLP`
- `DC`
- `CC`

The workflow target is `workflow/classical_baseline_benchmark.smk`.

## Species and Label Policy

Included species:

- `human`
- `celegans`
- `scerevisiae`
- `fgraminearum`

Label policy:

- `human`, `celegans`, and `scerevisiae` reuse the current ProGATE_v2 Phase 2A legacy replay label setup through `src.data.build_epgat_legacy_dataset.build_dataset`, which reads the legacy OGEE-backed labels for those species.
- `fgraminearum` is explicitly wired to the established new-label regime through:
  - `results/phase2b_new_label/labels/new_positive.tsv`
  - `results/phase2b_new_label/labels/new_negative.tsv`
- `old440` is not used anywhere in this workflow.

## Feature Settings

Feature combinations are defined over the legacy EPGAT feature blocks already migrated into ProGATE_v2:

- `ORT`
- `EXP`
- `SUB`
- `ORT_EXP`
- `ORT_SUB`
- `EXP_SUB`
- `ORT_EXP_SUB`

Method-by-method settings:

- `MLP`: all 7 feature combinations above.
- `RF`, `SVM`, `NB`: `ORT_EXP_SUB` only.
- `N2V_MLP`: `N2V` only.
- `DC`, `CC`: `network` only.

Benchmark design constraints implemented here:

- Tabular baselines use only tabular node features.
- `N2V_MLP` uses only Node2Vec embeddings.
- `DC` and `CC` use only graph topology from the PPI graph.

## Reused EPGAT Source Files

The migration reused logic and method definitions from these original EPGAT files:

- `/home/jiehuang/software/fungi/EPGAT/runners/run_mlp.py`
- `/home/jiehuang/software/fungi/EPGAT/runners/run_rf.py`
- `/home/jiehuang/software/fungi/EPGAT/runners/run_svm.py`
- `/home/jiehuang/software/fungi/EPGAT/runners/run_bayes.py`
- `/home/jiehuang/software/fungi/EPGAT/runners/run_n2v_mlp.py`
- `/home/jiehuang/software/fungi/EPGAT/runners/run_dc.py`
- `/home/jiehuang/software/fungi/EPGAT/runners/run_cc.py`
- `/home/jiehuang/software/fungi/EPGAT/runners/tools.py`

Adaptations made for ProGATE_v2:

- Data loading was moved onto `src.data.build_epgat_legacy_dataset`.
- Output layout was changed to the ProGATE_v2 benchmark directory structure.
- Fusarium labels were switched to the current new-label regime.
- `DC` and `CC` were kept as pure topology baselines to avoid feature leakage.
- Node2Vec was implemented through the `torch_geometric` Node2Vec class already available in the `EPGAT` environment.

## Workflow Files

Main workflow files:

- `workflow/classical_baseline_benchmark.smk`
- `scripts/run_classical_baseline_benchmark.sh`
- `configs/classical_baseline_benchmark.yaml`

Core implementation files:

- `src/train/run_classical_baseline_single.py`
- `src/graph/run_node2vec_embedding.py`
- `src/network/run_network_heuristics.py`
- `src/eval/aggregate_classical_baselines.py`
- `src/classical_baselines/common.py`

## Output Layout

Trainable methods write per-run outputs under:

- `outputs/classical_baseline_benchmark/{species}/{method}/{feature_setting}/run_{run_id}/`

Examples:

- `outputs/classical_baseline_benchmark/human/MLP/ORT/run_0/`
- `outputs/classical_baseline_benchmark/fgraminearum/RF/ORT_EXP_SUB/run_1/`
- `outputs/classical_baseline_benchmark/scerevisiae/N2V_MLP/N2V/run_2/`

Deterministic heuristics write to:

- `outputs/classical_baseline_benchmark/{species}/{method}/network/`

Examples:

- `outputs/classical_baseline_benchmark/celegans/DC/network/`
- `outputs/classical_baseline_benchmark/fgraminearum/CC/network/`

Final merged summary:

- `results/classical_baseline_benchmark/final_summary.tsv`

## How To Run

Shell launcher:

```bash
cd /home/jiehuang/software/fungi/ProGATE_v2
bash scripts/run_classical_baseline_benchmark.sh 8
```

The launcher runs Snakemake twice:

1. the full benchmark matrix
2. the final merged summary target

Manual Snakemake invocation:

```bash
cd /home/jiehuang/software/fungi/ProGATE_v2
export XDG_CACHE_HOME=/tmp
export MPLCONFIGDIR=/home/jiehuang/software/fungi/ProGATE_v2/.mplconfig
mapfile -t TARGETS < <(python - <<'PY'
import yaml

cfg = yaml.safe_load(open("configs/classical_baseline_benchmark.yaml"))
output_root = cfg["paths"]["benchmark_output_root"]
new_positive = cfg["paths"]["new_label_positive_path"]
new_negative = cfg["paths"]["new_label_negative_path"]
species = ["human", "celegans", "scerevisiae", "fgraminearum"]
trainable = {
    "MLP": ["ORT", "EXP", "SUB", "ORT_EXP", "ORT_SUB", "EXP_SUB", "ORT_EXP_SUB"],
    "RF": ["ORT_EXP_SUB"],
    "SVM": ["ORT_EXP_SUB"],
    "NB": ["ORT_EXP_SUB"],
    "N2V_MLP": ["N2V"],
}
heuristics = ["DC", "CC"]
run_ids = [0, 1, 2]

print(new_positive)
print(new_negative)
for sp in species:
    for method, feature_settings in trainable.items():
        for feature in feature_settings:
            for run_id in run_ids:
                print(f"{output_root}/{sp}/{method}/{feature}/run_{run_id}/metrics.tsv")
            print(f"{output_root}/{sp}/{method}/{feature}/aggregated_metrics.tsv")
for sp in species:
    for method in heuristics:
        print(f"{output_root}/{sp}/{method}/network/metrics.tsv")
        print(f"{output_root}/{sp}/{method}/network/aggregated_metrics.tsv")
PY
)

/home/jiehuang/anaconda3/envs/snakemake-latest/bin/snakemake \
  --snakefile workflow/classical_baseline_benchmark.smk \
  --cores 1 \
  --rerun-incomplete \
  "${TARGETS[@]}"

/home/jiehuang/anaconda3/envs/snakemake-latest/bin/snakemake \
  --snakefile workflow/classical_baseline_benchmark.smk \
  --cores 1 \
  --rerun-incomplete \
  results/classical_baseline_benchmark/final_summary.tsv
```

## Environment

The workflow uses:

- `/home/jiehuang/anaconda3/envs/EPGAT/bin/python` for the Python benchmark code
- `/home/jiehuang/anaconda3/envs/snakemake-latest/bin/snakemake` for workflow execution
