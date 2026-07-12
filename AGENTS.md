# Project identity

- The sole current project name is **EvoGATE**.
- EvoGATE is an evolution-aware essential-gene prediction and prioritization research framework for plant-pathogenic fungi.
- ProGATE and ProGATE_v2 are historical names. Preserve them only when describing historical artifacts, paths, or provenance.
- EPGAT is a methodological predecessor and source of inherited implementations. It is not an EvoGATE alias or internal product name.
- Bingo is an external comparison method and a historical source for some processed labels. It is not an EvoGATE component.
- The current primary biological target is *Fusarium graminearum* PH-1.
- The credible current application is computational essential-gene prioritization. RNA target discovery, off-target filtering, and dsRNA design are Planned.

# Scientific invariants

- Evolution-aware supervision is the first-class design principle.
- Label quality and label meaning take priority over model complexity.
- Evolution-aware label reconstruction is the primary scientific contribution.
- GraphSAGE is the main implementation carrier, not the sole innovation.
- PPI, ORT, EXP, SUB, and ESM2 are complementary evidence or representation blocks.
- Lethality, virulence, pathogenicity, fitness, and essentiality are not interchangeable terms.
- A transferred ortholog label is derived evidence, not direct target-species experimental validation.
- Essential-gene prediction is upstream of RNA target discovery.
- A predicted or prioritized candidate must never be described as experimentally validated.
- Candidate prioritization currently has no wet-lab validation.
- RNA target discovery, off-target filtering, and dsRNA design must remain labeled Planned until code, data, and validation exist.
- Do not optimize wording, thresholds, or experiments toward a desired result.

# Repository map

| Path | Responsibility | Handling rule |
|---|---|---|
| `configs/` | Frozen protocol and experiment configurations | Read before any run; do not silently edit |
| `data/manifests/` | Dataset inventories and source strategies | Provenance records |
| `data/derived_labels/` | Evolutionary transfer artifacts | Core evidence; read-only by default |
| `data/interim/` | Protocol bridges and preliminary evidence | Provenance/intermediate; read-only by default |
| `data/processed/` | Materialized labels, PPI, ORT, EXP, SUB, ESM2 | Data contract; do not overwrite |
| `docs/` | Project knowledge base and historical audits | Keep English/Chinese facts synchronized |
| `scripts/` | Figure builders and current/historical wrappers | Inspect paths before use |
| `src/data/` | Data adapters, bridges, label and graph loading | Main data implementation layer |
| `src/features/` | ESM2 and feature utilities | Feature implementation layer |
| `src/models/` | Classical, graph, and fusion models | Model implementation layer |
| `src/train/` | Training and experiment runners | Large runs require approval |
| `src/eval/`, `src/analysis/`, `src/plot/` | Metrics, summaries, interpretation, and Figures | Preserve claim-to-artifact traceability |
| `workflow/` | Snakemake workflows | Dry-run or execution requires approval when it can trigger work |
| `results/frozen_protocol/` | Frozen labels and splits | Immutable unless explicitly versioning a new protocol |
| `results/Figure*` | Figure summaries and retained artifacts | Do not overwrite |
| `results/label_rebuild_experiments/` | Historical reconstruction experiments | Historical; not a mainline runtime input |

# Canonical entry points

| Purpose | Canonical entry | Current status |
|---|---|---|
| Fusarium label materialization | `workflow/fgraminearum_label_materialization.smk` | Partially reproducible |
| Freeze protocol manifests | `python -m src.data.freeze_unified_protocol --config configs/frozen_protocol.yaml` | Implemented; writes frozen artifacts |
| Frozen benchmark | `workflow/frozen_protocol_benchmark.smk` | Implemented; large job |
| Single benchmark task | `python -m src.train.run_frozen_protocol_model` | Implemented |
| Frozen data loading | `src/data/frozen_protocol_loader.py` | Validated core loader |
| ESM2 extraction | `python -m src.features.extract_esm2_pooled` | Implemented; large job |
| ESM2 cache workflow | `workflow/prepare_esm2_cache.smk` | Implemented; large job |
| Candidate prioritization | `python -m src.eval.build_figure5_candidate_prioritization` | Blocked by missing `outputs/` inputs |
| Figure workflows | Named files under `workflow/Figure*.smk` | Experiment-specific; inspect before use |

Shell wrappers that hard-code `/home/jiehuang/software/fungi/ProGATE_v2` are Historical / non-portable. Do not present them as recommended EvoGATE entry points. Do not batch-edit them without a separate migration audit and approval.

# Frozen protocol

- Protocol version: `frozen_protocol_v1`.
- Main Fusarium protocol: `fgraminearum_newlabel`.
- Historical comparison protocol: `fgraminearum_oldlabel`.
- Formal split: 70% train, 10% validation, 20% test.
- Split seed: `20260409`.
- Training seeds: `1029`, `1030`, `1031`, `1032`, `1033`.
- Current newlabel: 1,097 positives and 10,868 negatives.
- Standard trainable-model threshold: `0.5`.
- A tuned threshold may be selected only from validation predictions.
- Never tune a threshold, model, feature set, epoch, graph threshold, or hyperparameter against the test split.
- All formal performance comparisons must use the same split and a compatible feature/evaluation contract.
- Do not silently change labels, sample membership, identifiers, split, seed, graph threshold, graph contract, features, model-selection metric, decision threshold, or aggregation method.
- Prefer AUPRC, MCC, and AUROC for this imbalanced task. Accuracy is supplementary.
- A five-seed mean and standard deviation describes optimization variability on one frozen split; it is not independent biological replication.
- Any new protocol must use a new explicit version and a separate output directory. It must not mutate `frozen_protocol_v1` artifacts.

# Evidence and status language

Use these status labels exactly:

| Status | Meaning |
|---|---|
| Implemented | Code or workflow exists for the stated behavior |
| Validated | Existing evidence verifies the stated contract or result within its declared scope |
| Partially implemented | Some components exist, but the end-to-end capability is incomplete |
| Partially validated | Evidence exists but is incomplete, conflicting, or insufficient for the full claim |
| Planned | Intended future work without current implementation evidence |
| Unknown | Available evidence cannot establish the fact |
| Blocked | A required dependency, source, environment, or artifact is missing |
| Historical | Retained for provenance, replay, or migration, not current mainline use |

- Every important statement should point to a config, source module, manifest, result table, or report.
- Do not infer behavior from a file name alone.
- Do not convert a goal, Figure title, old report narrative, or comment into an implemented result.
- Where multiple result versions conflict, state that multiple conflicting artifacts exist and link `docs/INCONSISTENCIES.md`.
- Do not select the largest or most favorable result.
- Use “predicted candidate” or “prioritized candidate,” not “validated target.”
- Feature-group zero-out perturbation is not SHAP, GNNExplainer, or causal attribution.

# Working rules

- Audit the affected flow before modifying it: config, caller, producer, consumer, and output contract.
- Search for existing helpers and patterns before adding code.
- Keep changes scoped to the requested module and scientific contract.
- Preserve raw evidence, processed data, frozen manifests, and published result artifacts.
- New runs must use an independent, versioned output directory.
- Never overwrite a frozen result or reuse its directory for a changed protocol.
- Record the exact command, working directory, environment, resolved config, input paths, output path, protocol version, split version, seed, and exit status for approved runs.
- Preserve source identifiers and canonical identifiers together in audit outputs.
- Fail visibly on unresolved alignment where silent loss could change the sample universe.
- Do not silently fall back from newlabel to oldlabel or from real ESM2 to mock embeddings.
- Keep historical replay entry points available until their replacement has been validated.
- Performance comparisons require identical splits and clearly reported feature/model contracts.
- Do not change test-set interpretation after seeing results.
- Update English and Chinese project documentation together when scientific facts change.

# Validation requirements

For documentation-only changes:

- verify every referenced repository path exists, unless explicitly labeled missing or Planned
- compare English and Chinese headings, key numbers, paths, status labels, tables, and Mermaid nodes
- check Markdown heading hierarchy, fenced blocks, tables, relative links, and Mermaid accessibility fields
- report line counts, byte sizes, and SHA-256 checksums when Git metadata is unavailable

For small code changes, after explicit authorization:

- run the smallest import, `--help`, parser, or focused self-check that exercises the change
- verify label counts, split counts, feature dimensions, graph node/edge counts, and output schema if affected
- confirm the change does not modify the frozen split or sample universe
- write test outputs only to a new temporary or versioned directory

For protocol, data, model, or Figure changes:

- obtain user confirmation before execution
- state expected compute and overwrite behavior
- run a dry-run or lightweight preflight first when available
- compare against the authoritative frozen contract
- retain resolved configs, logs, predictions, metrics, and checksums
- document any deviation before interpreting performance

# Naming and migration

- EvoGATE is the only current project name.
- ProGATE and ProGATE_v2 are historical names.
- EPGAT remains the predecessor name; do not rename inherited references when they identify provenance.
- Bingo remains an external project name.
- Historical absolute paths may be evidence. Report them before changing them.
- Do not perform global search-and-replace on project names or paths.
- Do not rewrite generated historical reports only to modernize names.
- A path migration must be separated from scientific or algorithmic changes.
- Stop migration work if sample counts, labels, mappings, splits, features, edges, metrics, or outputs change.
- Follow `docs/MIGRATION_GUIDE.md` and record old-to-new mappings explicitly.

# Safety

- Treat `data/` and `results/` as read-only unless the user explicitly authorizes writes.
- Do not run full training, benchmark, Snakemake, ESM2 extraction, graph rebuild, Figure regeneration, or large analysis without explicit approval.
- Do not run `rm`, `git reset`, bulk `mv`, or overwriting `cp` operations.
- Do not delete, move, rename, or replace historical artifacts during routine work.
- Do not create symlinks in the repository root.
- Do not download data or models or install dependencies without approval.
- Do not load large model, embedding, matrix, or binary artifacts for routine inspection.
- Inspect large `data/` and `results/` trees by directory, metadata, headers, and small samples.
- Verify output destinations before an approved run; reject existing non-versioned destinations.
- Never expose credentials, tokens, private endpoints, or personal environment secrets in code or documentation.

# Documentation standards

- Main research-facing documents must have English and Chinese versions with matching structure, numbers, paths, status labels, table rows, and diagrams.
- Code comments should primarily be in English.
- Audit reports may be Chinese-first when requested.
- Use repository-relative paths for current instructions.
- Preserve absolute paths only when documenting migration or provenance.
- Use Mermaid for architecture and flow diagrams; include `accTitle` and `accDescr` where supported.
- Distinguish scientific contributions, engineering contributions, limitations, and future work.
- Keep authoritative facts in one clear location and link instead of duplicating unstable detail.
- Put operational commands in the runbook, scientific logic in scientific design, data schemas in the data dictionary, and current blockers in project status/inconsistencies.
- Do not use promotional claims such as “state-of-the-art,” “universally applicable,” “experimentally validated RNA targets,” or “ready-to-use dsRNA platform” without direct evidence.

# Known blockers

- `.git/` is empty, so repository history and normal Git status/diff checks are unavailable.
- The repository-local `outputs/` directory is absent; some summaries and Figure5 inputs cannot be rebuilt from this workspace.
- No environment lock, dependency manifest, container, or release installation procedure exists.
- Several runtime and data-builder paths are tied to historical Linux or macOS locations.
- Some evaluation bytecode exists without corresponding source modules.
- The upstream producer and exact rule for `weak_positive_confidence` are missing.
- A complete standalone exclusion manifest is absent.
- Figure2a and Figure2b contain conflicting values for nominally similar settings.
- The old label-scarcity narrative conflicts with the available ranking table.
- Gated fusion improves some metrics and reduces others; it is not uniformly superior.
- Standard-label source descriptions conflict between OGEE-oriented narratives and Bingo/EPGAT manifest provenance.
- Old protocol audit documents cite some entry points that are missing or renamed.

Read `docs/INCONSISTENCIES.md`, `docs/PROJECT_STATUS.md`, and `docs/REPRODUCIBILITY.md` before making or evaluating scientific claims.

