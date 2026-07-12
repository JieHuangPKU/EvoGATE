# EvoGATE project status

_Evidence-based implementation and validation matrix._

---

## Status vocabulary

This document uses only: **Implemented**, **Validated**, **Partially implemented**, **Partially validated**, **Planned**, **Unknown**, **Blocked**, and **Historical**.

## Module status matrix

| Module | Scientific purpose | Implementation status | Validation status | Canonical entry point | Expected output | Current blocker | Evidence path |
|---|---|---|---|---|---|---|---|
| Yeast essentiality transfer artifact | Supply evolutionary supervision | Partially implemented | Partially validated | Unknown upstream generator | `ph1_yeast_essential_ortholog_labels.tsv` | Confidence generator missing | `data/derived_labels/` |
| PH-1 canonical bridge | Resolve protein and legacy IDs | Implemented | Validated | `src.data.build_fgraminearum_newlabel_bridge` | Bridge and unresolved audits | External configured sources affect rerun | `data/processed/essential_gene/fgraminearum/bridge/` |
| Lethal evidence preparation | Select direct lethal support | Implemented | Validated | `src.data.prepare_fgraminearum_label_materialization_sources` | `lethal_positive_gene_list.tsv` | Some configured inputs are historical paths | `data/interim/protocol_refactor/fgraminearum_label_materialization/` |
| Newlabel materialization | Define current supervision | Implemented | Validated | `src.data.materialize_fgraminearum_label_regimes` | 1,097 positives; 10,868 negatives | Full upstream reconstruction incomplete | `data/processed/essential_gene/fgraminearum/newlabel/` |
| Oldlabel replay | Preserve historical comparison | Implemented | Historical | Same materialization module | oldlabel tables | Not a mainline regime | `data/processed/essential_gene/fgraminearum/oldlabel/` |
| Frozen labels and splits | Prevent model-specific resampling | Implemented | Validated | `src.data.freeze_unified_protocol` | Frozen label and split manifests | Writes existing frozen results | `results/frozen_protocol/` |
| Orthology features | Represent evolutionary context | Implemented | Partially validated | `src.data.build_inparanoid_ortholog_matrix` | `orthologs.csv` | Builder has historical machine paths | `data/processed/OR/` |
| Expression features | Represent expression profiles | Implemented | Partially validated | `src.data.build_expression_profile_csv` | `profile.csv` | Raw/source build portability incomplete | `data/processed/EXP/` |
| Localization features | Represent cellular context | Implemented | Partially validated | `src.data.build_subloc_csv_from_compartments` | `subloc.csv` | Raw/source build portability incomplete | `data/processed/LC/` |
| ESM2 embeddings | Represent protein sequence | Implemented | Validated | `src.features.extract_esm2_pooled` | `esm2_pooled.pt` | Local model and environment paths | `data/processed/ESM2/` |
| STRING graph | Define default PPI topology | Implemented | Validated | `src.data.frozen_protocol_loader` | Filtered edge table/index | Complete raw rebuild not portable | `data/processed/PPI/` |
| eFG graph comparison | Test graph-source robustness | Implemented | Partially validated | `workflow/Figure4_graph_robustness.smk` | Figure4 summaries | Large rerun requires approval | `results/Figure4/` |
| Unified frozen loader | Align labels, graph, and features | Implemented | Validated | `src.data.frozen_protocol_loader` | In-memory benchmark bundle | Environment not locked | `src/data/frozen_protocol_loader.py` |
| Classical baselines | Provide non-graph comparisons | Implemented | Partially validated | `src.train.run_frozen_protocol_model` | Per-run metrics/predictions | `outputs/` missing | `results/Figure1/summary/` |
| Topology baselines | Measure network-only information | Implemented | Partially validated | Same runner | node2vec/DC/CC metrics | Backend/environment reproducibility | `results/Figure1/summary/` |
| GNN families | Compare message-passing models | Implemented | Partially validated | Same runner | GAT/GCN/GIN/GraphSAGE runs | Conflicting Figure2 artifacts | `results/Figure2a/`, `results/Figure2b/` |
| ESM2 comparison | Quantify sequence-representation contribution | Implemented | Partially validated | `workflow/Figure3a_fusarium_graphsage_esm2_comparison.smk` | Figure3a summaries | Per-run outputs absent | `results/Figure3a/` |
| ESM2 dimension ablation | Compare truncated embedding dimensions | Implemented | Partially validated | `workflow/Figure3b_fusarium_graphsage_esm2_dim_ablation.smk` | Figure3b summaries | Wrapper non-portable | `results/Figure3b/` |
| Fusion ablation | Compare concatenation and gates | Implemented | Partially validated | Figure3c workflows | Figure3c summaries | Mixed metric direction; claims unresolved | `results/Figure3c*/` |
| Label-scarcity benchmark | Test reduced-label behavior | Implemented | Partially validated | `workflow/label_scarcity_benchmark.smk` | Scarcity metrics and plots | Old narrative conflicts with ranking | `results/Figure2_label_scarcity/` |
| Graph robustness | Test PPI thresholds and sources | Implemented | Partially validated | `workflow/Figure4_graph_robustness.smk` | Figure4 tables and plots | Statistical claim scope incomplete | `results/Figure4/` |
| Representation analysis | Examine hidden and input spaces | Implemented | Partially validated | `workflow/Figure5_representation_mechanism.smk` | UMAP and summary artifacts | Depends on missing run outputs for rebuild | `results/Figure5/` |
| Group zero-out analysis | Estimate feature-group sensitivity | Implemented | Partially validated | `workflow/Figure5d_feature_group_attribution.smk` | Group perturbation summaries | Not causal attribution or SHAP | `results/Figure5/` |
| Candidate prioritization | Rank experimental candidates | Partially implemented | Partially validated | `src.eval.build_figure5_candidate_prioritization` | Candidate rank tables | Missing `outputs/`; no wet-lab validation | `results/Figure5_new_candidate_prioritization/` |
| RNA target discovery | Evaluate RNA target suitability | Planned | Unknown | None | None | No implementation | None |
| Off-target filtering | Remove host/non-target matches | Planned | Unknown | None | None | No implementation | None |
| dsRNA design | Design silencing sequences | Planned | Unknown | None | None | No implementation | None |

## Repository maturity

EvoGATE is a research repository with a validated supervision core and substantial experimental artifacts. Release-grade portability and reproducibility are **Blocked**. It must not be described as a production-ready software release or validated RNA-target platform.
