# EvoGATE 数据字典

_当前标签、特征、图、评价和候选优选流程使用的可验证核心表与 artifact。_

---

## 约定

只列出能够从现有 header、code schema 或 config 验证的 column。标为 **Unknown** 的值不会被推断。历史文本 artifact 中可能同时出现 empty string 与 `NA`；consumer 应遵循 producer module，不应假定全局 missing-value rule。

## 核心 artifact 字典

| Artifact | Path pattern | 目的 | Row entity | Key identifier | Major columns | Missing-value convention | Provenance | Producer | Consumer | Status |
|---|---|---|---|---|---|---|---|---|---|---|
| Standard labels | `data/processed/essential_gene/<species>/labels.standard.tsv` | 模式生物 benchmark label 和历史 Fusarium label | Gene | `gene_id` 或 file-specific canonical field | Label 和 source/audit fields；完整 schema 因文件而异 | Empty/`NA`；见文件 | Dataset manifest 和 source audit | Historical dataset audit/build code | Frozen protocol | Historical / Validated artifact |
| Newlabel labels | `data/processed/essential_gene/fgraminearum/newlabel/labels.tsv` | 当前 Fusarium supervision | Canonical gene | `canonical_gene_id` | `graph_gene_id`、`label`、`label_text`、`regime`、source fields | TSV read 后 empty string | Lethal list、transfer table、bridge、evidence mirror | `materialize_fgraminearum_label_regimes.py` | Frozen protocol | Validated |
| Positive labels | `.../newlabel/positive_genes.tsv` | Positive-set provenance | Canonical gene | `canonical_gene_id` | `graph_gene_id`、`label`、`positive_sources`、`construction_bucket`、support fields | Empty string | PHI lethal 与 high-confidence transfer | 同一 materializer | Label freeze、audit | Validated |
| Negative labels | `.../newlabel/negative_genes.tsv` | Negative-set provenance | Canonical gene | `canonical_gene_id` | `graph_gene_id`、`label`、`construction_bucket`、`source_manifest` | Empty string | 经 exclusion 的 resolved none pool | 同一 materializer | Label freeze、audit | Validated |
| Materialized split | `.../{newlabel,oldlabel}/split.tsv` | Label-regime split | Labeled gene | `canonical_gene_id` / `graph_gene_id` | `label`、`split`、`split_seed`、`split_strategy`、`split_version` | 不应缺失 split | Materialized label | 同一 materializer | Frozen protocol | Validated |
| Frozen labels | `results/frozen_protocol/labels/*.tsv` | Model-independent label contract | Labeled gene | `canonical_gene_id`、`graph_gene_id` | Label、regime、provenance、protocol fields | Empty string | Processed label source | `freeze_unified_protocol.py` | Frozen loader | Validated |
| Frozen splits | `results/frozen_protocol/splits/*.tsv` | Shared split contract | Labeled gene | `graph_gene_id` | `canonical_gene_id`、`label`、`split`、`split_seed`、`split_version` | 不应缺失 split | Frozen labels | `freeze_unified_protocol.py` | Frozen loader | Validated |
| Protein bridge | `data/processed/essential_gene/fgraminearum/bridge/protein_to_canonical_bridge.tsv` | 映射 PH-1 protein 到 canonical gene | Source protein | `source_protein_id` | `resolved_canonical_gene_id`、`bridge_status`、`bridge_method`；support columns | Unresolved target 为空 | Unified map、legacy map、sequence/header evidence | `build_fgraminearum_newlabel_bridge.py` | Label materializer | Validated |
| Source mapping | `.../bridge/source_to_canonical_mapping.tsv` | 合并 evidence 中的 source ID | Source identifier | `source_id` | `canonical_gene_id`；mapping evidence columns | Empty string | Bridge input 与 modality map | 同一 bridge builder | Source preparation | Validated |
| High transfer candidates | `.../bridge/high_confidence_yeast_transfer_candidates.tsv` | Resolved transfer-positive component | Transfer-supported gene/protein | `canonical_gene_id` 与 `ph1_gene_id` | `orthogroup_id`、`bridge_method`、confidence support | Empty string | Yeast transfer 与 bridge | 同一 bridge builder | Label materializer | Validated |
| Yeast transfer | `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv` | Evolutionary essentiality evidence | PH-1 protein | `ph1_gene_id` | `orthogroup_id`、yeast support flags/IDs、occupancy/copy fields、`weak_positive_confidence` | Support ID 中出现 `NA` | Summary 记录的 OrthoFinder result 与 yeast list | Unknown upstream generator | Bridge 与 materializer | Partially implemented |
| ORT features | `data/processed/OR/<species>/orthologs.csv` | Orthology feature matrix | Gene | Source 中 `Gene`，标准化为 `legacy_gene_id` | Target-species binary columns | Loader 将 numeric missing 填零 | InParanoid/orthoXML build report | `build_inparanoid_ortholog_matrix.py` | Frozen loader | Partially reproducible |
| EXP features | `data/processed/EXP/<species>/profile.csv` | Expression feature matrix | Gene | Source 中 `Gene`，标准化为 `legacy_gene_id` | Expression/sample columns | Processed/loader path 将 missing 填零 | GEO/GTEx build report | `build_expression_profile_csv.py` 与 update script | Frozen loader | Partially reproducible |
| SUB features | `data/processed/LC/<species>/subloc.csv` | Localization matrix | Gene | Source 中 `Gene`，标准化为 `legacy_gene_id` | Controlled localization categories | Binary absent/unknown 记为零 | COMPARTMENTS/eFG build report | `build_subloc_csv_from_compartments.py` | Frozen loader | Partially reproducible |
| ESM2 embeddings | `data/processed/ESM2/<species>/esm2_pooled.pt` | Protein sequence representation | Protein/gene key | Embedding dictionary key | `embeddings`、optional `metadata`；当前 full vector dimension 为 1,280 | Frozen loader 中 missing key 为 fatal | Species protein FASTA 和 local ESM2 model | `extract_esm2_pooled.py` | Frozen loader | Validated artifact |
| STRING PPI | `data/processed/PPI/<species>/string.csv` | 默认 graph source | Interaction edge | Composite `A`、`B` | `A`、`B`、可用时为 `combined_score` | 移除 invalid/empty endpoint | STRING v12 mapping/build report | `build_string_csv_from_string_v12.py` | Frozen loader | Validated artifact |
| eFG PPI | `data/processed/PPI/fgraminearum/eFG_ppis.txt` | Fusarium graph-source comparison | Interaction edge | Source-specific endpoints | Confidence 与 endpoint fields；精确 schema 见 adapter | Unknown | eFG | `prepare_figure4_efg_graph.py` 进行适配 | Figure4 workflow | Partially validated |
| Per-run predictions | `outputs/<experiment>/<protocol>/<model>/<feature>/run_<seed>/predictions.tsv` | Node-level score | Graph node | `graph_gene_id` | `canonical_gene_id`、`split`、`label`、`pred_score`、`pred_label`、protocol/model fields | Unlabeled node 的 label 为空 | Frozen bundle 与 trained model | `run_frozen_protocol_model.py` | Evaluation/Figure/candidate module | Blocked in workspace |
| Per-run metrics | `outputs/.../metrics.tsv` | Run-level validation/test metric 与 provenance | One run | protocol/model/feature/seed | `val_*`、`test_*`、count、contract、path、feature dimension | 不适用字段为空 | Prediction 与 run config | 同一 model runner | Aggregator | Blocked in workspace |
| Aggregated metrics | `results/Figure*/**/*aggregated_metrics.tsv` | 跨 seed 描述性 summary | Model/feature/protocol group | Composite fields | Mean、standard deviation、run count、seed list | 不可用 metric 可能为 `NA` | Per-run metric | Figure-specific aggregator | Figure/manuscript | Partially validated |
| Candidate ranking | `results/Figure5_new_candidate_prioritization/Figure5_new_candidate_rank_table.tsv` | 比较 baseline/ESM2 ranking | Graph gene | `gene_id` | Score、rank、percentile、delta、label fields、top-k flags | Unlabeled `true_label` 可缺失 | 跨 seed Figure3a prediction | `build_figure5_candidate_prioritization.py` | Candidate review | Partially implemented |

## ID 规则

- Fusarium canonical ID：`fgraminearum::FGRAMPH1_*`
- Fusarium graph-facing ID：通常为 `FGRAMPH1_*`
- Protein transfer source：`XP_*`
- Orthogroup：`OG*`
- 模式生物 ID 因 source 而异，由 processed feature builder 标准化

ID conversion 必须显式且可审计。外观相似的 gene、transcript、protein 与 orthogroup ID 不可互换。

## 状态与 schema 缺口

历史 standard-label table 和 raw eFG input 的完整 schema 因 artifact 而异，尚未形成单一 contract，因此本文不进行补全。`weak_positive_confidence` 的 upstream algorithm 为 **Unknown**。缺失的 `outputs/` tree 阻止直接验证所有当前 per-run table instance。

