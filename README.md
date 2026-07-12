# EvoGATE

_Evolution-aware essential-gene prediction and prioritization for plant-pathogenic fungi._

---

> EvoGATE is currently a research repository, not yet a production-ready software release.

## Project summary

EvoGATE reconstructs auditable essential-gene supervision from phenotype and evolutionary evidence, combines that supervision with protein interaction networks, biological features, and ESM2 embeddings, and prioritizes candidate essential genes in *Fusarium graminearum*.

## Current scientific scope

The current scope is essential-gene prediction and computational candidate prioritization. The primary target is *F. graminearum* PH-1. Human, *Caenorhabditis elegans*, *Saccharomyces cerevisiae*, and *Drosophila melanogaster* are retained as benchmark species.

### What EvoGATE currently does

- **Validated**: materializes an evolution-aware `fgraminearum_newlabel` regime with 1,097 positives and 10,868 negatives
- **Validated**: freezes a stratified 70/10/20 train/validation/test split
- **Implemented**: integrates PPI, orthology (ORT), expression (EXP), subcellular localization (SUB), and ESM2 representations
- **Implemented**: compares MLP, RF, SVM, NB, node2vec, GAT, GCN, GIN, GraphSAGE, and network heuristics
- **Partially implemented**: ranks genome-wide candidates and identifies ESM2-associated rank changes

### What EvoGATE does not yet do

- **Planned**: RNA target discovery
- **Planned**: host and non-target off-target filtering
- **Planned**: dsRNA or siRNA design
- **Blocked**: release-grade end-to-end reproduction in this workspace

Candidate predictions are computational hypotheses, not experimentally validated essential genes or effective RNA targets.

## Core scientific contribution

The primary contribution is evolution-aware label reconstruction. PHI-supported lethal evidence is combined with high-confidence transfer from yeast essential orthologs after explicit PH-1 identifier bridging. GraphSAGE is the main modeling vehicle, not the sole innovation.

## Repository status warning

The repository was migrated from the historical names ProGATE and ProGATE_v2. EPGAT is a methodological predecessor and source of inherited implementations; Bingo is an external comparison project and a historical source for some standard labels. Neither is an EvoGATE alias or internal component.

Several wrappers still contain historical absolute paths, the repository-local `outputs/` directory is absent, `.git/` is empty, environment lock files are absent, and some source or upstream generation steps are missing. See [Known inconsistencies](docs/INCONSISTENCIES.md) and [Reproducibility](docs/REPRODUCIBILITY.md).

## Canonical workflow entry points

| Purpose | Canonical entry | Status |
|---|---|---|
| Label materialization | `workflow/fgraminearum_label_materialization.smk` | Partially reproducible |
| Frozen benchmark | `workflow/frozen_protocol_benchmark.smk` | Implemented; large job |
| Single benchmark task | `python -m src.train.run_frozen_protocol_model` | Implemented |
| ESM2 cache preparation | `workflow/prepare_esm2_cache.smk` | Implemented; large job |
| Candidate prioritization | `python -m src.eval.build_figure5_candidate_prioritization` | Blocked by missing `outputs/` inputs |

Historical Shell wrappers that hard-code `/home/jiehuang/software/fungi/ProGATE_v2` are non-portable and are not recommended entry points.

## Frozen protocol

| Setting | Value |
|---|---|
| Protocol | `frozen_protocol_v1` |
| Split | 70% train / 10% validation / 20% test |
| Split seed | `20260409` |
| Training seeds | `1029`, `1030`, `1031`, `1032`, `1033` |
| Main metrics | AUPRC, MCC, AUROC |
| Threshold policy | Fixed `0.5` for standard runs; any tuned threshold must be selected on validation data |

## Repository map

| Path | Role |
|---|---|
| `configs/` | Frozen protocol and experiment configuration |
| `data/` | Manifests, derived labels, interim provenance, and processed modalities |
| `docs/` | Project knowledge base, migration audits, and scientific documentation |
| `scripts/` | Thin and historical execution wrappers plus Figure builders |
| `src/` | Data, feature, graph, model, training, evaluation, and analysis modules |
| `workflow/` | Snakemake workflows |
| `results/` | Frozen manifests, summaries, Figures, and historical experiment artifacts |

## Documentation index

- [Project overview](docs/PROJECT_OVERVIEW.md)
- [Scientific design](docs/SCIENTIFIC_DESIGN.md)
- [Architecture](docs/ARCHITECTURE.md)
- [Data flow](docs/DATA_FLOW.md)
- [Runbook](docs/RUNBOOK.md)
- [Project status](docs/PROJECT_STATUS.md)
- [Reproducibility](docs/REPRODUCIBILITY.md)
- [Data dictionary](docs/DATA_DICTIONARY.md)
- [Project history](docs/PROJECT_HISTORY.md)
- [Migration guide](docs/MIGRATION_GUIDE.md)
- [Manuscript mapping](docs/MANUSCRIPT_MAPPING.md)
- [Glossary](docs/GLOSSARY.md)
- [Known inconsistencies](docs/INCONSISTENCIES.md)
- [Chinese README](README.zh-CN.md)

## Reproducibility, citation, and license status

- **Reproducibility**: Partially reproducible; not release-grade
- **Citation**: No canonical citation file has been established
- **License**: No repository license has been established
- **Project maturity**: Research repository with validated core artifacts and unresolved portability and provenance blockers

