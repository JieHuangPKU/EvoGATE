# EvoGATE known inconsistencies

_Conflicts and blockers that must remain visible until resolved by evidence and explicit review._

---

## Interpretation policy

An inconsistency is not resolved by choosing the newest, largest, or most favorable number. Resolution requires a traceable input, code version, configuration, run artifact, and documented decision. Until then, affected claims are **Partially validated**, **Unknown**, or **Blocked**.

## Issue register

| ID | Issue | Affected files | Observed evidence | Authoritative interpretation | Impact | Required resolution | Status |
|---|---|---|---|---|---|---|---|
| INC-001 | 80/10/10 versus 70/10/20 split | `REVIEW_REPORT.md`; `configs/frozen_protocol.yaml`; `src/data/freeze_unified_protocol.py`; frozen split summaries | Review report states 80/10/10; config, code, split version, and counts encode 70/10/20 | The frozen protocol is 70/10/20 | Incorrect methods reporting and sample counts | 已修复 — REVIEW_REPORT.md corrected to 70/10/10 on 2026-07-12 | 已修复 |
| INC-002 | Figure2a and Figure2b report different values for nominally similar settings | `results/Figure2a/`; `results/Figure2b/`; corresponding configs/workflows | GraphSAGE `ORT_EXP_SUB` summaries differ across artifacts | Multiple conflicting artifacts exist; timestamp evidence confirms different code versions (Fig2a: Apr 12, Fig2b: Apr 11, model runner: May 1) | Feature/model claims may depend on run lineage | Recover per-run outputs, resolved configs, code versions, and compare contracts | Partially validated |
| INC-003 | Label-scarcity narrative conflicts with actual ranking | `REVIEW_REPORT.md`; `results/Figure2_label_scarcity/summary/label_scarcity_report.md`; ranking table | Older text claims GraphSAGE superiority; at 10% labels MLP ranks first and GraphSAGE third by AUPRC | The available ranking does not support broad GraphSAGE superiority | Invalid headline claim if repeated | 部分修复 — REVIEW_REPORT.md Table 0.2 已加注实际排名 (2026-07-12)；label_scarcity_report.md 仍需审计 | Partially validated |
| INC-004 | Gated fusion improves some metrics and reduces others | `results/Figure3c*/Figure3c*_final_summary.md`; old narrative | Gated AUPRC increases relative to concatenation in one comparison, while AUROC, MCC, and F1 decrease; residual variants are not uniformly superior | Fusion effects are metric-dependent | Cannot claim gated fusion is universally better | Predefine primary endpoint, use paired tests, and report trade-offs | Partially validated |
| INC-005 | OGEE versus Bingo/EPGAT label-source descriptions | `REVIEW_REPORT.md`; `data/manifests/essential_gene_dataset_manifest.tsv`; migration docs | Some descriptions call standard labels OGEE-based; manifest records Bingo or Bingo+EPGAT strategies | Frozen tables are authoritative runtime inputs; upstream source attribution remains mixed | Provenance and external-credit risk | 部分修复 — REVIEW_REPORT.md "OGEE标准"已改为"Bingo标准" (2026-07-12)；上游Bingo→OGEE lineage仍待确认 | Partially validated |
| INC-006 | EvoGATE identity conflicts with ProGATE_v2 naming | package docstrings, configs, scripts, docs, result reports | Current formal name is EvoGATE; internal text and paths still use ProGATE_v2 | EvoGATE is the only current name; ProGATE/ProGATE_v2 are Historical | User confusion and non-portable entry points | Perform reviewed, path-by-path migration in a later phase | Historical |
| INC-007 | Repository-local `outputs/` is missing | Figure5 summaries, workflows, candidate builder | Results cite `outputs/Figure3a/...`, but `outputs/` does not exist in this workspace | Existing summaries are retained artifacts, not fully rebuildable local evidence | Breaks result-to-run traceability | Recover or regenerate outputs in a versioned directory after approval | Blocked |
| INC-008 | `.git/` is empty | `.git/`; failed Git repository inspection | Directory exists but contains no repository metadata | No local Git history is available | Cannot verify authorship, versions, or diffs | 已修复 — git init, baseline commit 3649532, pushed to GitHub on 2026-07-12 | 已修复 |
| INC-009 | Environment lock is absent | Repository root and configs | No `environment.yml`, lock file, `requirements.txt`, or `pyproject.toml` exists; configs contain machine paths | Environment is not release-reproducible | Imports and numerical behavior may vary | Capture tested environments and lock dependencies in a separately approved phase | Blocked |
| INC-010 | Bytecode exists without source | `src/eval/__pycache__/evaluate_graph_model*.pyc`, `evaluate_baseline*.pyc`, `evaluate_support_graph_baseline*.pyc` | `.pyc` artifacts are present while corresponding `.py` source is absent | Missing source cannot be treated as auditable implementation | Evaluation lineage may be incomplete | Recover source from history or regenerate results through auditable modules | Blocked |
| INC-011 | Yeast-transfer confidence generator is missing | `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`; summary; consuming modules | Artifact contains `weak_positive_confidence`; repository has no verified producer for this field | Downstream consumption is known; confidence assignment is Unknown | Core evolutionary rule cannot be fully independently rebuilt | Recover producer, formalize rules, and create provenance/checksum record | Blocked |
| INC-012 | Old documents cite missing or renamed entry points | `docs/protocol_refactor/protocol_state_audit.md` and migration documents | Documents name workflows/scripts/evaluation modules absent from current inventory | These documents are Historical audits, not current runbooks | Misleading developer guidance | Map each historical entry to current, archived, or missing status | Historical |

## Additional material inconsistencies

| ID | Issue | Interpretation | Status |
|---|---|---|---|
| INC-013 | Generated reports embed `/home/jiehuang/software/fungi/ProGATE_v2` output paths | Preserve as historical provenance; do not treat as current root | Historical |
| INC-014 | Data builders contain `/Users/jiehuang/.../ProGATE_v2` paths | Processed artifacts exist, but raw rebuild portability is blocked | Blocked |
| INC-015 | `data/processed/ESM2-bk/` duplicates current ESM2 caches | Do not delete or select without checksum and consumer audit | 已修复 — ESM2-bk directory does not exist in current workspace (verified 2026-07-12) | 已修复 |
| INC-016 | Existing Figure claims often report mean and standard deviation without paired inference | Descriptive evidence exists; inferential claim strength is limited | Partially validated |

## Resolution order

1. Establish authoritative protocol and result-release identifiers.
2. Recover Git history, missing source, environment definition, and run-level outputs.
3. Reconstruct label-source and confidence provenance.
4. Recompute only after explicit approval, into a new versioned output directory.
5. Update manuscript claims after paired, contract-matched comparisons.

