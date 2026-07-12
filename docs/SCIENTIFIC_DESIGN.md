# EvoGATE scientific design

_Scientific hypotheses, supervision design, representations, evaluation, interpretation, and limitations._

---

## Scientific question

Can essential genes be predicted and prioritized in a label-scarce plant-pathogenic fungus by reconstructing supervision from lethal phenotype evidence and conserved essential orthology, then integrating that supervision with interaction networks, biological features, and protein-language representations?

## Core hypothesis

EvoGATE tests three linked hypotheses:

1. Label reconstruction changes the learnable biological problem more substantially than replacing one GNN family with another.
2. PPI topology, evolutionary features, expression, localization, and ESM2 contain complementary information.
3. Stable multi-seed predictions can prioritize experimental candidates, while remaining computational hypotheses until independently validated.

## Evolution-aware label reconstruction

### Positive definition

The `fgraminearum_newlabel` positive set is:

```text
PHI-supported lethal canonical genes
UNION
resolved high-confidence yeast-essentiality-transfer canonical genes
```

The materialized counts are 77 lethal genes, 1,045 transfer-supported genes, 25 overlapping genes, and 1,097 unique positives. Evidence is consumed by `src/data/materialize_fgraminearum_label_regimes.py` under `configs/fgraminearum_label_materialization.yaml`.

### Negative definition

The negative set begins with successfully bridged genes whose transfer confidence is `none`. Genes are removed if they have virulence or pathogenicity evidence or occur in the positive set. The materialized negative set contains 10,868 genes.

### Exclusion logic

The following do not become supervised negatives:

- high-confidence positives
- genes with virulence or pathogenicity evidence
- low- or medium-confidence transfer genes
- unresolved or ambiguous identifier mappings
- other genes not admitted by the positive or negative definitions

The logic is present across source and materialization code, but a complete standalone exclusion table is absent. Status: **Partially implemented**.

### Canonical ID bridge

`src/data/build_fgraminearum_newlabel_bridge.py` maps PH-1 protein-level identifiers into canonical `fgraminearum::FGRAMPH1_*` gene space. It uses the configured unified map, legacy mapping evidence, sequence/header evidence, and modality-specific mapping tables. Unresolved high-confidence rows are retained in an audit output rather than silently admitted.

### Evolutionary evidence

`data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv` records *S. cerevisiae* and *S. pombe* essential-ortholog support, fungal occupancy, copy statistics, single-copy indicators, support classes, and confidence. The script that originally assigned `weak_positive_confidence` is missing; its generation is **Unknown**, while the downstream artifact and its use are documented.

## Feature modalities

| Modality | Input | Representation | Role |
|---|---|---|---|
| ORT | `data/processed/OR/<species>/orthologs.csv` | Binary orthology matrix | Evolutionary context |
| EXP | `data/processed/EXP/<species>/profile.csv` | Numeric expression profile | Condition-associated activity |
| SUB | `data/processed/LC/<species>/subloc.csv` | Binary localization matrix | Cellular context |
| ESM2 | `data/processed/ESM2/<species>/esm2_pooled.pt` | Mean-pooled protein embedding | Sequence-derived representation |
| Degree | PPI edge table | Scalar graph feature | Local topology |

Feature normalization is computed from training nodes in `src/data/frozen_protocol_loader.py`. ESM2 alignment is strict: missing graph-node embeddings raise an error rather than causing silent row loss.

## Graph construction

The default graph source is `data/processed/PPI/<species>/string.csv`. STRING edges are filtered at `combined_score >= 300` under the frozen protocol, self-loops are removed, duplicates are dropped, and the default graph contract is `undirected_symmetrized`. Edge weights are disabled in the frozen configuration. Figure4 additionally evaluates STRING thresholds and eFG graph sources.

## Model families

| Family | Models | Scientific role |
|---|---|---|
| Tabular | MLP, RF, SVM, NB | Non-graph feature baselines |
| Topology | node2vec+MLP, degree, closeness | Network-only or embedding baselines |
| Graph neural network | GAT, GCN, GIN, GraphSAGE | Message-passing comparisons |
| Fusion variants | Concatenation, gated, residual gated, weighted BCE variants | Multimodal design ablation |

GraphSAGE is the primary implementation carrier. Existing artifacts do not justify treating it as the sole scientific contribution.

## Frozen evaluation

| Item | Frozen value |
|---|---|
| Split | 70% train / 10% validation / 20% test |
| Split seed | `20260409` |
| Training seeds | `1029`-`1033` |
| Standard decision threshold | `0.5` |
| Tuned-threshold source | Validation split only |
| Primary metrics | AUPRC, MCC, AUROC |

AUPRC is emphasized because the positive class is much smaller than the negative class. MCC summarizes all four confusion-matrix cells and remains informative under imbalance. AUROC is reported as a complementary ranking metric. Test data must not be used for model, threshold, or hyperparameter selection.

## Ablation and interpretation

Implemented analyses include oldlabel versus newlabel comparison, feature combinations, ESM2 inclusion, ESM2 dimension truncation, fusion variants, label scarcity, graph thresholds and sources, hidden-representation analysis, and feature-group zero-out perturbation.

Feature-group attribution is not SHAP or GNNExplainer. It measures the change after setting selected standardized columns to zero in a frozen model. This supports group-level sensitivity, not causal residue-level interpretation.

## Candidate ranking

`src/eval/build_figure5_candidate_prioritization.py` combines baseline and ESM2 predictions, ranks genes, summarizes top-k overlap and rank changes, and profiles candidate groups. Current candidate tables mix labeled and unlabeled graph nodes where indicated by their label fields; downstream users must filter explicitly for their intended experimental question.

Status: **Partially implemented** and not wet-lab validated.

## Contributions and boundaries

### Scientific contributions

- Evolution-aware supervision for fungal essentiality
- Auditable positive, negative, and exclusion rules
- Frozen comparison of complementary biological representations
- Candidate prioritization grounded in prediction stability and evidence profiles

### Engineering contributions

- Canonical identifier bridge
- Shared frozen loader and model contract
- Per-run provenance, prediction, metric, and feature-schema outputs
- Structured result aggregation and Figure workflows

### Current limitations

- Missing upstream confidence generator
- No explicit complete exclusion manifest
- Conflicting result versions for some comparisons
- No release-grade statistical inference across all claims
- Missing run-level `outputs/` in the current workspace
- No experimental validation of prioritized candidates

### Future work

RNA target discovery, off-target filtering, conserved target-region selection, dsRNA design, cross-species transfer evaluation, and experimental validation are **Planned**. They are not current EvoGATE capabilities.

