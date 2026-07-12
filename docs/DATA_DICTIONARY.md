# EvoGATE data dictionary

_Verified core tables and artifacts used by the current label, feature, graph, evaluation, and prioritization flows._

---

## Conventions

Only columns confirmed from existing headers, code schemas, or configuration are listed. A value marked **Unknown** is intentionally not inferred. Empty strings and `NA` may both occur in historical text artifacts; consumers must follow the producing module rather than assume a universal missing-value rule.

## Core artifact dictionary

| Artifact | Path pattern | Purpose | Row entity | Key identifier | Major columns | Missing-value convention | Provenance | Producer | Consumer | Status |
|---|---|---|---|---|---|---|---|---|---|---|
| Standard labels | `data/processed/essential_gene/<species>/labels.standard.tsv` | Benchmark labels for model organisms and historical Fusarium | Gene | `gene_id` or file-specific canonical field | Label and source/audit fields; exact full schema varies | Empty/`NA`; see file | Dataset manifest and source audits | Historical dataset audit/build code | Frozen protocol | Historical / Validated artifact |
| Newlabel labels | `data/processed/essential_gene/fgraminearum/newlabel/labels.tsv` | Current Fusarium supervision | Canonical gene | `canonical_gene_id` | `graph_gene_id`, `label`, `label_text`, `regime`, source fields | Empty string after TSV reads | Lethal list, transfer table, bridge, evidence mirror | `materialize_fgraminearum_label_regimes.py` | Frozen protocol | Validated |
| Positive labels | `.../newlabel/positive_genes.tsv` | Positive-set provenance | Canonical gene | `canonical_gene_id` | `graph_gene_id`, `label`, `positive_sources`, `construction_bucket`, support fields | Empty string | PHI lethal and high-confidence transfer | Same materializer | Label freeze, audits | Validated |
| Negative labels | `.../newlabel/negative_genes.tsv` | Negative-set provenance | Canonical gene | `canonical_gene_id` | `graph_gene_id`, `label`, `construction_bucket`, `source_manifest` | Empty string | Resolved none pool after exclusions | Same materializer | Label freeze, audits | Validated |
| Materialized split | `.../{newlabel,oldlabel}/split.tsv` | Label-regime split | Labeled gene | `canonical_gene_id` / `graph_gene_id` | `label`, `split`, `split_seed`, `split_strategy`, `split_version` | No missing split expected | Materialized labels | Same materializer | Frozen protocol | Validated |
| Frozen labels | `results/frozen_protocol/labels/*.tsv` | Model-independent label contract | Labeled gene | `canonical_gene_id`, `graph_gene_id` | Label, regime, provenance, protocol fields | Empty string | Processed label sources | `freeze_unified_protocol.py` | Frozen loader | Validated |
| Frozen splits | `results/frozen_protocol/splits/*.tsv` | Shared split contract | Labeled gene | `graph_gene_id` | `canonical_gene_id`, `label`, `split`, `split_seed`, `split_version` | No missing split expected | Frozen labels | `freeze_unified_protocol.py` | Frozen loader | Validated |
| Protein bridge | `data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv` | Map PH-1 proteins to canonical genes | Source protein | `source_protein_id` | `resolved_canonical_gene_id`, `bridge_status`, `bridge_method`; support columns | Empty string for unresolved targets | Unified map, legacy maps, sequence/header evidence | `build_fgraminearum_newlabel_bridge.py` | Label materializer | Validated |
| Source mapping | `.../bridge/source_to_canonical_mapping.tsv` | Consolidate source IDs across evidence | Source identifier | `source_id` | `canonical_gene_id`; mapping evidence columns | Empty string | Bridge inputs and modality maps | Same bridge builder | Source preparation | Validated |
| High transfer candidates | `.../bridge/high_confidence_yeast_transfer_candidates.tsv` | Resolved transfer-positive component | Transfer-supported gene/protein | `canonical_gene_id` and `ph1_gene_id` | `orthogroup_id`, `bridge_method`, confidence support | Empty string | Yeast transfer plus bridge | Same bridge builder | Label materializer | Validated |
| Yeast transfer | `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv` | Evolutionary essentiality evidence | PH-1 protein | `ph1_gene_id` | `orthogroup_id`, yeast support flags/IDs, occupancy/copy fields, `weak_positive_confidence` | `NA` appears in support IDs | OrthoFinder results and yeast lists recorded in summary | Unknown upstream generator | Bridge and materializer | Partially implemented |
| ORT features | `data/processed/OR/<species>/orthologs.csv` | Orthology feature matrix | Gene | `Gene` in source, normalized to `legacy_gene_id` | Target-species binary columns | Numeric missing values filled with zero by loader | InParanoid/orthoXML build reports | `build_inparanoid_ortholog_matrix.py` | Frozen loader | Partially reproducible |
| EXP features | `data/processed/EXP/<species>/profile.csv` | Expression feature matrix | Gene | `Gene` in source, normalized to `legacy_gene_id` | Expression/sample columns | Missing values filled with zero in processed/loader paths | GEO/GTEx build reports | `build_expression_profile_csv.py` and update scripts | Frozen loader | Partially reproducible |
| SUB features | `data/processed/LC/<species>/subloc.csv` | Localization matrix | Gene | `Gene` in source, normalized to `legacy_gene_id` | Controlled localization categories | Binary absent/unknown represented as zero in matrix | COMPARTMENTS/eFG build reports | `build_subloc_csv_from_compartments.py` | Frozen loader | Partially reproducible |
| ESM2 embeddings | `data/processed/ESM2/<species>/esm2_pooled.pt` | Protein sequence representation | Protein/gene key | Embedding dictionary key | `embeddings`, optional `metadata`; vector dimension 1,280 in current full cache | Missing keys are fatal in frozen loader | Species protein FASTA and local ESM2 model | `extract_esm2_pooled.py` | Frozen loader | Validated artifact |
| STRING PPI | `data/processed/PPI/<species>/string.csv` | Default graph source | Interaction edge | Composite `A`, `B` | `A`, `B`, `combined_score` when available | Invalid/empty endpoints removed | STRING v12 mapping/build reports | `build_string_csv_from_string_v12.py` | Frozen loader | Validated artifact |
| eFG PPI | `data/processed/PPI/fgraminearum/eFG_ppis.txt` | Fusarium graph-source comparison | Interaction edge | Source-specific endpoints | Confidence and endpoint fields; exact schema in adapter | Unknown | eFG | `prepare_figure4_efg_graph.py` adapts it | Figure4 workflow | Partially validated |
| Per-run predictions | `outputs/<experiment>/<protocol>/<model>/<feature>/run_<seed>/predictions.tsv` | Node-level scores | Graph node | `graph_gene_id` | `canonical_gene_id`, `split`, `label`, `pred_score`, `pred_label`, protocol/model fields | Unlabeled nodes have empty label | Frozen bundle and trained model | `run_frozen_protocol_model.py` | Evaluation/Figure/candidate modules | Blocked in workspace |
| Per-run metrics | `outputs/.../metrics.tsv` | Run-level validation/test metrics and provenance | One run | protocol/model/feature/seed | `val_*`, `test_*`, counts, contract, paths, feature dimension | Empty for inapplicable fields | Predictions and run config | Same model runner | Aggregators | Blocked in workspace |
| Aggregated metrics | `results/Figure*/**/*aggregated_metrics.tsv` | Across-seed descriptive summaries | Model/feature/protocol group | Composite fields | Means, standard deviations, run counts, seed list | `NA` may appear for unavailable metrics | Per-run metrics | Figure-specific aggregators | Figures/manuscript | Partially validated |
| Candidate ranking | `results/Figure5_new_candidate_prioritization/Figure5_new_candidate_rank_table.tsv` | Compare baseline and ESM2 rankings | Graph gene | `gene_id` | Scores, ranks, percentiles, deltas, label fields and top-k flags | Unlabeled `true_label` may be missing | Figure3a predictions across seeds | `build_figure5_candidate_prioritization.py` | Candidate review | Partially implemented |

## Identifier rules

- Canonical Fusarium ID: `fgraminearum::FGRAMPH1_*`
- Graph-facing Fusarium ID: usually `FGRAMPH1_*`
- Protein transfer source: `XP_*`
- Orthogroup: `OG*`
- Model-organism identifiers vary by source and are normalized by processed feature builders

Identifier conversion must be explicit and audited. Similar-looking gene, transcript, protein, and orthogroup identifiers are not interchangeable.

## Status and schema gaps

The full schema of historical standard-label tables and raw eFG input is not normalized in this dictionary because it varies by artifact and has not been established as a single contract. The upstream algorithm for `weak_positive_confidence` is **Unknown**. The absent `outputs/` tree prevents direct verification of all current per-run table instances.

