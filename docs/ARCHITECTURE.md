# EvoGATE architecture

_Code architecture, scientific architecture, module boundaries, and execution relationships._

---

## Seven-layer code architecture

```mermaid
flowchart TD
    accTitle: EvoGATE code layers
    accDescr: Seven repository layers convert evidence into labels, features, graphs, trained models, evaluations, and candidate-ranking artifacts.

    knowledge_layer["Knowledge layer<br/>data/derived_labels, data/interim, data/manifests"]
    label_layer["Label layer<br/>data/processed/essential_gene, results/frozen_protocol"]
    feature_layer["Feature layer<br/>data/processed/OR, EXP, LC, ESM2"]
    graph_layer["Graph layer<br/>data/processed/PPI, src/data, src/graph"]
    model_layer["Model layer<br/>src/models, src/train"]
    evaluation_layer["Evaluation layer<br/>src/eval, src/analysis, src/plot"]
    application_layer["Application layer<br/>Figure5 candidate prioritization"]

    knowledge_layer --> label_layer
    knowledge_layer --> feature_layer
    label_layer --> graph_layer
    feature_layer --> graph_layer
    graph_layer --> model_layer --> evaluation_layer --> application_layer
```

| Layer | Responsibility | Main paths | Canonical entry |
|---|---|---|---|
| Knowledge | Evidence, species scope, transfer artifacts, ID provenance | `data/derived_labels/`, `data/interim/`, `data/manifests/` | `src.data.build_fgraminearum_newlabel_bridge` |
| Label | Positive/negative regimes and frozen splits | `data/processed/essential_gene/`, `results/frozen_protocol/` | `workflow/fgraminearum_label_materialization.smk`, `src.data.freeze_unified_protocol` |
| Feature | ORT, EXP, SUB, ESM2 feature blocks | `data/processed/OR/`, `data/processed/EXP/`, `data/processed/LC/`, `data/processed/ESM2/`, `src/features/` | Modality-specific builders; no unified feature workflow |
| Graph | PPI filtering, node universe, edge index, topology embeddings | `data/processed/PPI/`, `src/data/frozen_protocol_loader.py`, `src/graph/` | `src.data.frozen_protocol_loader` |
| Model | Classical, topology, GNN, and fusion models | `src/models/`, `src/train/` | `src.train.run_frozen_protocol_model` |
| Evaluation | Metrics, aggregation, ablation, interpretation, plots | `src/eval/`, `src/analysis/`, `src/plot/`, `workflow/` | Figure workflows and evaluation modules |
| Application | Candidate ranking and future target discovery | `src/eval/build_figure5_candidate_prioritization.py`, `results/Figure5*` | Candidate module; RNA layer is Planned |

## Scientific architecture

```mermaid
flowchart TD
    accTitle: EvoGATE scientific architecture
    accDescr: Phenotype and molecular evidence are transformed through evolutionary supervision and graph learning into candidate rankings, with RNA and dsRNA stages explicitly marked as future work.

    evidence["Evidence<br/>PHI phenotypes, yeast essentiality, molecular data"]
    evolution["Evolution<br/>orthogroups, conservation, copy structure"]
    labels["Labels<br/>positive, negative, exclusion, frozen split"]
    representation["Representation<br/>ORT, EXP, SUB, ESM2"]
    graph_learning["Graph learning<br/>STRING/eFG and model families"]
    evaluation["Evaluation<br/>frozen metrics, ablation, robustness"]
    prioritization["Candidate prioritization<br/>scores, stability, rank shifts"]
    rna_target["RNA target discovery<br/>Planned"]
    dsrna_design["dsRNA design<br/>Planned"]

    evidence --> evolution --> labels
    evidence --> representation
    labels --> graph_learning
    representation --> graph_learning --> evaluation --> prioritization
    prioritization -. future .-> rna_target -. future .-> dsrna_design
```

| Scientific stage | Status |
|---|---|
| Evidence assembly | Partially implemented |
| Evolutionary transfer artifact | Partially implemented |
| Label materialization | Validated |
| Multimodal representation | Validated |
| Graph learning | Validated |
| Evaluation | Partially validated |
| Candidate prioritization | Partially implemented |
| RNA target discovery | Planned |
| dsRNA design | Planned |

## Main execution chain

```mermaid
flowchart LR
    accTitle: Frozen benchmark call chain
    accDescr: The frozen protocol workflow materializes manifests, loads aligned graph data, runs one model task per configuration and seed, then aggregates metrics and produces plots.

    frozen_config["configs/frozen_protocol.yaml"]
    frozen_workflow["workflow/frozen_protocol_benchmark.smk"]
    freeze_protocol["src.data.freeze_unified_protocol"]
    frozen_loader["src.data.frozen_protocol_loader"]
    model_runner["src.train.run_frozen_protocol_model"]
    aggregate_runs["src.eval.aggregate_frozen_protocol_runs"]
    plot_rules["workflow/plots.smk"]

    frozen_config --> frozen_workflow
    frozen_workflow --> freeze_protocol
    frozen_workflow --> model_runner
    model_runner --> frozen_loader
    frozen_workflow --> aggregate_runs --> plot_rules
```

## Configuration model

`configs/frozen_protocol.yaml` defines repository-relative data roots, protocol names, frozen runtime settings, feature roots, ESM2 caches, label sources, model families, and model hyperparameters. Figure-specific YAML files reference this base configuration and override experiment scope or model variants.

Resolved runtime configurations are expected in per-run output directories. Current `results/Figure3*/runtime/` artifacts preserve some rendered configs, but the primary `outputs/` tree is missing from this workspace.

## Data and identifier contracts

The primary Fusarium canonical identifier is `fgraminearum::FGRAMPH1_*`; graph-facing files may use the prefix-stripped `FGRAMPH1_*` value. `frozen_protocol_loader.py` joins labels, graph nodes, feature rows, and ESM2 keys through explicit graph and canonical identifiers.

The loader constructs a node universe from the union of graph and labeled nodes, maps the frozen split to node indices, normalizes numeric features using training nodes, and returns a single bundle consumed by model families.

## Model and output contracts

A standard run writes `predictions.tsv`, `metrics.tsv`, `feature_schema.tsv`, `edge_table.tsv`, `split_manifest.tsv`, `resolved_config.yaml`, and model-specific artifacts such as `best_model.pt`, `model.pkl`, `training_log.tsv`, or ESM2 alignment audits.

The current workspace contains aggregated results and Figures but lacks the main `outputs/` tree. Consequently, this output contract is **Implemented** in code but **Blocked** for complete local reconstruction.

## Legacy boundaries

- `docs/epgat_migration/` records EPGAT migration history
- `docs/protocol_refactor/` records ProGATE_v2 protocol refactoring
- several `scripts/run_*.sh` files hard-code the historical ProGATE_v2 path
- legacy training and data adapters remain under `src/` for controlled replay
- historical artifacts must not be treated as current canonical entry points unless explicitly named in this document

See [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) for the non-destructive migration policy.

## Dependencies and portability

The code imports Python, pandas, NumPy, PyYAML, scikit-learn, PyTorch, graph libraries, Snakemake, R, and plotting packages. No authoritative environment lock exists. Some configured Python and cache paths are machine-specific. Dependency and hardware reproducibility are therefore **Blocked** at release grade.
