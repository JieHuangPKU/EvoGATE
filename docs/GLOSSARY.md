# EvoGATE glossary

_Controlled terminology for scientific, data, model, and application documentation._

---

## Biological terms

| Term | Definition in EvoGATE |
|---|---|
| Essential gene | A gene required for viability under a stated biological context. EvoGATE labels operationalize this concept from defined evidence; a prediction is not proof of essentiality. |
| Lethal phenotype | A mutant or perturbation phenotype indicating non-viability under the recorded assay. It is the direct phenotype category used for the PHI-supported positive component. |
| Virulence | The degree of damage or disease severity caused by a pathogen. Reduced virulence does not by itself imply loss of viability. |
| Pathogenicity | The ability to cause disease. Loss of pathogenicity is not equivalent to lethality. |
| Fitness | Relative reproductive or growth performance in a defined environment. Reduced fitness is broader than essentiality. |
| Core essential | A gene required across a broad range of conditions or conserved as essential across lineages. EvoGATE does not experimentally establish universal core essentiality. |
| Context-specific essentiality | Essentiality restricted to a condition, host, tissue, developmental stage, or environmental context. |
| RNA target | A transcript or gene selected for potential RNA-mediated suppression after essentiality, accessibility, conservation, and off-target assessment. EvoGATE has not implemented this stage. |
| dsRNA design | Selection of a double-stranded RNA sequence intended to generate effective silencing fragments while controlling off-target risk. Status: Planned. |

## Evolution and supervision

| Term | Definition in EvoGATE |
|---|---|
| Evolution-aware supervision | Supervision constructed using explicit evolutionary relationships and phenotype evidence rather than treating sparse target-species annotations as the only labels. |
| Ortholog transfer | Transfer of evidence from a gene in one species to an orthologous target-species gene under a stated confidence rule. It is derived evidence, not direct target-species validation. |
| Single-copy ortholog | An orthologous relationship or orthogroup with one relevant copy per species under the specified scope. The exact scope must be reported. |
| Orthogroup | A group of genes descended from a common ancestral gene, represented by identifiers such as `OG*` in the transfer artifact. |
| Positive label | A gene admitted as essential under the declared positive rule: PHI-supported lethal evidence or resolved high-confidence yeast-essentiality transfer. |
| Negative label | A resolved `none`-confidence gene retained after virulence/pathogenicity and positive-set exclusions. It is an operational negative, not universal proof of non-essentiality. |
| Exclusion set | Genes not admitted as positives or negatives because of uncertain evidence, biological exclusions, unresolved mapping, or other protocol rules. A complete standalone manifest is currently absent. |
| Weak positive confidence | Existing categorical field in the yeast-transfer artifact. Its downstream use is known; its upstream producer is missing and therefore Unknown. |

## Protocol and evaluation

| Term | Definition in EvoGATE |
|---|---|
| Frozen protocol | Versioned labels, split, graph/feature contracts, seeds, metrics, and settings that must remain fixed across formal comparisons. |
| Split seed | `20260409`, used to materialize the stratified 70/10/20 split. |
| Training seed | One of `1029`-`1033`, used for model stochasticity on the frozen split. |
| AUPRC | Area under the precision-recall curve. It emphasizes positive-class ranking and is primary for the imbalanced essential-gene task. |
| AUROC | Area under the receiver operating characteristic curve. It measures ranking across true-positive and false-positive rates and is complementary to AUPRC. |
| MCC | Matthews correlation coefficient, a threshold-dependent summary using all four confusion-matrix cells. |
| Validation threshold | A decision threshold selected only from validation predictions, never from test outcomes. |
| Fixed threshold | The standard trainable-model decision threshold of `0.5` in the frozen runner. |

## Representation and model terms

| Term | Definition in EvoGATE |
|---|---|
| PPI | Protein-protein interaction network used as graph structure. Default processed source: STRING; eFG is used in source-robustness analysis. |
| ORT | Orthology feature block derived from processed cross-species orthology matrices. |
| EXP | Numeric expression feature block. |
| SUB | Binary subcellular localization feature block. |
| ESM2 embedding | Mean-pooled protein representation generated from the configured ESM2 model; the current full embedding dimension is 1,280. |
| GraphSAGE | The primary message-passing implementation carrier in current experiments. It is not the sole EvoGATE innovation. |
| Gated fusion | A learned fusion variant for omics and ESM2 blocks. Existing results show metric-dependent trade-offs. |
| Group zero-out perturbation | Inference-time masking of standardized feature columns to estimate group sensitivity. It is not SHAP, GNNExplainer, or causal attribution. |

## Application and status terms

| Term | Definition in EvoGATE |
|---|---|
| Candidate prioritization | Computational ranking using prediction scores, cross-seed behavior, rank changes, and evidence profiles. It does not imply experimental validation. |
| Predicted candidate | A gene prioritized by a computational model. Use this term instead of “validated target.” |
| Implemented | Code or workflow exists for the stated behavior. |
| Validated | Existing evidence verifies the stated contract or result within its declared scope. |
| Partially implemented | Some required components exist, but the end-to-end capability is incomplete. |
| Partially validated | Evidence exists but is incomplete, conflicting, or insufficient for the full claim. |
| Planned | Intended future work without current implementation evidence. |
| Unknown | Available evidence cannot establish the fact. |
| Blocked | A required dependency or artifact is missing. |
| Historical | Retained for provenance, replay, or migration rather than current mainline use. |

