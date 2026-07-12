# EvoGATE manuscript mapping

_Traceability from scientific claims to evidence, code, configuration, data, outputs, and Figures._

---

## Claim-to-artifact matrix

| Claim | Evidence | Code | Config | Input data | Output artifact | Figure | Status | Caveat |
|---|---|---|---|---|---|---|---|---|
| Evolution-aware label reconstruction produces the current Fusarium supervision regime | Materialized counts and source manifests | `build_fgraminearum_newlabel_bridge.py`; `prepare_fgraminearum_label_materialization_sources.py`; `materialize_fgraminearum_label_regimes.py` | `fgraminearum_label_materialization.yaml` | Transfer table, evidence mirror, unified/legacy maps | `data/processed/essential_gene/fgraminearum/newlabel/` | Figure1B/1C | Validated | Upstream confidence producer missing |
| Newlabel contains 1,097 positives and 10,868 negatives | `summary.tsv`, frozen label summary | Same materializer and freeze module | Label and frozen configs | Materialized positive/negative tables | Newlabel and frozen summaries | Figure1C | Validated | Exclusion set is not a standalone complete manifest |
| Newlabel and oldlabel differ in composition and model behavior | Old/new comparison tables | Figure2a aggregation and feature-combo runner | `Figure2a_fusarium_label_compare_graphsage.yaml` | Oldlabel/newlabel frozen manifests | `results/Figure2a/` | Figure2a | Partially validated | Figure2 artifact versions conflict; see `INCONSISTENCIES.md` |
| Multiple classical, topology, and GNN families were benchmarked | Frozen model list and Figure1 summary | `run_frozen_protocol_model.py`; aggregator | `frozen_protocol.yaml` | Frozen labels/splits and processed features | `results/Figure1/summary/` | Figure1 benchmark | Partially validated | Per-run `outputs/` missing; no complete paired inference |
| Biological feature groups contribute differently | Feature-combination results | Feature-combo runner and Figure2b workflow | `Figure2b_fusarium_gnn_feature_ablation.yaml` | ORT, EXP, SUB, graph | `results/Figure2b/` | Figure2b | Partially validated | Conflicts with Figure2a for nominally similar settings |
| ESM2 adds information to GraphSAGE in Fusarium | Four-setting comparison summary | Figure3a workflow and aggregation | `Figure3a_fusarium_graphsage_esm2_comparison.yaml` | ESM2 cache, ORT, EXP, SUB, graph | `results/Figure3a/` | Figure3a | Partially validated | Descriptive five-seed evidence; per-run outputs absent |
| ESM2 dimensional truncation changes performance | Dimension summaries | Figure3b workflow | `Figure3b_fusarium_graphsage_esm2_dim_ablation.yaml` | Full ESM2 cache and frozen bundle | `results/Figure3b/` | Figure3b | Partially validated | Truncates embedding dimensions; does not compare ESM2 layers/models |
| Gated fusion changes performance relative to concatenation | Figure3c summaries | Fusion model and Figure3c workflows | Figure3c YAML files | ORT/EXP/SUB/ESM2 frozen bundle | `results/Figure3c*/` | Figure3c | Partially validated | AUPRC and threshold metrics move in different directions; no universal superiority claim |
| Performance under reduced labels was measured | Ranking and per-run scarcity tables | Split builder, workflow, summarizer | Base frozen config through scarcity workflow | Frozen newlabel split and features | `results/Figure2_label_scarcity/` | Label-scarcity Figure | Partially validated | MLP ranks first at 10%; old GraphSAGE-superiority narrative is contradicted |
| PPI threshold and source affect performance | Threshold and source comparison tables | Figure4 graph adapter, runner, summarizer | `frozen_protocol.yaml` Figure4 section | STRING/eFG graphs and frozen bundle | `results/Figure4/` | Figure4 | Partially validated | Primarily Fusarium; inference scope limited |
| ESM2 changes hidden representations and feature-group sensitivity | UMAP, hidden summaries, zero-out tables | `src/analysis/` Figure5 modules | Figure5 workflows plus Figure3 runtime configs | Predictions, checkpoints, feature schemas | `results/Figure5/` | Figure5a-5d | Partially validated | Zero-out analysis is not SHAP/GNNExplainer or causal attribution |
| EvoGATE prioritizes candidate essential genes | Candidate ranks, top-k overlap, rank-shift tables | `build_figure5_candidate_prioritization.py` | CLI/runtime arguments | Figure3a predictions across seeds | `results/Figure5_new_candidate_prioritization/` | Figure5 | Partially implemented | Missing upstream `outputs/`; candidates lack wet-lab validation |
| EvoGATE discovers effective RNA targets | None | None | None | None | None | None | Planned | This claim must not be made |

## Recommended paper structure

| Section | Main evidence |
|---|---|
| Background and problem | Figure1A and label-source summaries |
| Scientific hypothesis | Frozen protocol and planned claim hierarchy |
| Evolution-aware label reconstruction | Figure1B/1C and label materialization artifacts |
| Framework | Loader, representations, graph, and model architecture |
| Benchmark | Figure1, Figure2, and Figure4 with conflict disclosure |
| ESM2 | Figure3 and Figure5 group perturbation |
| Interpretability | Figure5 representation and perturbation limitations |
| Candidate prioritization | Figure5 ranks as computational hypotheses |
| Discussion | Label provenance, imbalance, graph scope, portability, and missing validation |
| Future RNA target discovery | Planned only |

## Claim language policy

Use “implemented” for code paths, “validated” only where existing evidence supports the stated contract, “predicted” or “prioritized” for candidates, and “planned” for RNA applications. Do not use “experimentally validated” for candidate rankings. Claims affected by conflicting artifacts must cite [INCONSISTENCIES.md](INCONSISTENCIES.md) rather than selecting a favorable version.

