# EvoGATE project overview

_Scientific scope, contribution, application boundary, and current maturity._

---

## Biological problem

Experimental essential-gene evidence is sparse in plant-pathogenic fungi. Lethality, virulence, pathogenicity, and reduced fitness are related but non-equivalent phenotypes. Identifier changes between proteins, transcripts, legacy annotations, and canonical PH-1 genes further complicate label construction.

EvoGATE addresses the resulting supervision problem before applying machine learning. Its current target is *Fusarium graminearum* PH-1, with four model organisms retained for comparative benchmarks.

## Label scarcity

The mainline Fusarium regime contains 1,097 positive and 10,868 negative labels. Although this is larger than the historical oldlabel regime, positives remain a minority. Accuracy alone can therefore conceal poor positive-class performance; AUPRC and MCC are primary metrics.

The existing label-scarcity experiment is **Validated as an artifact** but its older narrative is inconsistent with the actual ranking. At 10% retained labels, MLP ranks first by AUPRC and GraphSAGE ranks third. No broad GraphSAGE-superiority claim should be made from this experiment.

## Evolution-aware supervision

The current positive set is the union of two evidence components:

1. 77 PHI-supported lethal genes that satisfy the protocolized evidence and canonical-ID requirements
2. 1,045 high-confidence yeast-essentiality-transfer genes after PH-1 protein-to-gene bridging

Twenty-five genes occur in both components, producing 1,097 unique positives. Negatives are resolved genes with `weak_positive_confidence == none`, after excluding genes with virulence or pathogenicity evidence and all positives. The final negative count is 10,868.

Low- and medium-confidence transfers, unresolved mappings, and biologically excluded genes do not enter the supervised label set. A single explicit exclusion manifest is not currently available, so the exclusion layer is **Partially implemented**.

## Multimodal graph learning

EvoGATE combines complementary inputs:

| Modality | Biological role | Current status |
|---|---|---|
| PPI | Graph structure for neighborhood aggregation | Implemented |
| ORT | Cross-species orthology presence and evolutionary context | Implemented |
| EXP | Expression profiles | Implemented |
| SUB | Subcellular localization indicators | Implemented |
| ESM2 | Protein-sequence representation | Implemented |

GraphSAGE is the main model carrier. GAT, GCN, GIN, MLP, RF, SVM, NB, node2vec, degree centrality, and closeness centrality provide comparison points. The scientific claim concerns supervision and complementary evidence, not the uniqueness of GraphSAGE.

## Candidate prioritization

Genome-wide predictions, cross-seed rank summaries, ESM2-associated rank changes, and feature-group perturbation analyses have been computed. These artifacts support **Partially implemented** candidate prioritization.

The candidates have not been validated by gene deletion, conditional knockdown, infection assays, or RNA interference experiments. They must be described as predicted or prioritized candidates.

## Intended application

The near-term application is prioritizing essential-gene hypotheses for experimental follow-up in *F. graminearum*. RNA target discovery, off-target filtering, and dsRNA design are **Planned** downstream applications. No current module establishes RNA target suitability or dsRNA efficacy.

## Contributions

### Scientific contributions

- Evolution-aware reconstruction of a fungal essential-gene label regime
- Explicit separation of lethal evidence, transferred evidence, negatives, and exclusions
- Evaluation of network, omics, evolutionary, and protein-language representations under a frozen protocol
- Computational prioritization of candidate essential genes

### Engineering contributions

- Auditable PH-1 identifier bridge and source manifests
- Frozen label and split contracts shared across model families
- Unified multimodal loader and structured per-run outputs
- Reusable Figure aggregation and interpretation workflows

## Current status

| Area | Status | Interpretation |
|---|---|---|
| Label reconstruction | Validated | Counts, sources, and materialized tables exist |
| Frozen evaluation contract | Validated | Split and seeds are explicit and materialized |
| Model and feature workflows | Implemented | Code and result summaries exist |
| Main performance claims | Partially validated | Conflicting result versions and incomplete statistical testing remain |
| Candidate prioritization | Partially implemented | Computed but not experimentally validated |
| Release-grade reproduction | Blocked | Missing outputs, locks, source modules, and portable entry points |
| RNA target discovery | Planned | No implementation evidence |

## Limitations

The upstream generator that assigned `weak_positive_confidence` is missing. Some evaluation modules exist only as bytecode. The current workspace lacks `outputs/`, and several workflows refer to historical absolute paths. Multiple result artifacts disagree for nominally similar Figure2 settings. These limitations are tracked in [INCONSISTENCIES.md](INCONSISTENCIES.md).

