# EvoGATE Snakemake Style Guide

_Maintenance standard for `workflow/*.smk` files._

---

## 1. Rule Naming

### Convention: `verb_object` in lowercase snake_case

| Good | Avoid | Reason |
|------|-------|--------|
| `run_trainable_model` | `materialize_or_run_figure2a` | Compound verbs with "or" obscure intent |
| `build_label_scarcity_splits` | `figure5a` | Noun-only names don't describe action |
| `summarize_figure4` | `run_figure3cc_old_gated_wbce` | Abbreviated suffixes (cc, ca) hard to parse |
| `prepare_efg_graph` | `Figure2a` | Capital letters break convention |
| `extract_esm2_pooled_embeddings` | — | Verb_noun pattern, clear and searchable |
| `plot_label_scarcity` | — | Standard plot prefix |
| `aggregate_frozen_protocol` | — | Standard aggregate prefix |

### Exception

`rule all:` is the standard Snakemake default target name. Keep as-is.

### Current violations (workflow/*.smk)

| File | Rule | Issue | Suggested rename |
|------|------|-------|-----------------|
| Figure2a | `materialize_or_run_figure2a` | Compound verb, vague | `run_figure2a_graphsage_comparison` |
| Figure2b | `materialize_or_run_figure2b` | Same | `run_figure2b_gnn_feature_ablation` |
| Figure3c | `materialize_or_run_figure3c` | Same | `run_figure3c_gated_fusion` |
| Figure3cA | `materialize_or_run_figure3ca` | Same + abbreviation | `run_figure3ca_gated_residual` |
| Figure3cB | `materialize_or_run_figure3cb` | Same + abbreviation | `run_figure3cb_gated_residual_wbce` |
| Figure5 | `figure5a` ... `figure5d` | Noun-only | `build_figure5a_umap_error_transition` etc. |

---

## 2. Configuration and Path Management

### Rule: No hardcoded paths in workflow files

**Current violations:**

```python
# workflow/fgraminearum_label_materialization.smk:4
PYTHON_BIN = "/home/jiehuang/anaconda3/bin/python"  # ❌ HARDCODED

# workflow/label_scarcity_benchmark.smk:197
~/anaconda3/bin/Rscript src/plot/plot_label_scarcity.R  # ❌ HARDCODED
```

**Required fix:**

```python
# Move to configs/fgraminearum_label_materialization.yaml
runtime:
  python_bin: "/home/jiehuang/anaconda3/bin/python"

# Move to configs/frozen_protocol.yaml
runtime:
  rscript_bin: "~/anaconda3/bin/Rscript"
```

### Rule: All paths derived from configs or Path variables

```python
# Good
OUTPUT_ROOT = config["paths"]["figure1_benchmark_output_root"]
SUMMARY_DIR = f"{RESULTS_ROOT}/summary"

# Acceptable (legacy — migrate to f-strings over time)
LABEL_SUMMARY = FIGURE1_LABELS_DIR + "/label_protocol_summary.md"

# Bad
PYTHON_BIN = "/home/jiehuang/anaconda3/bin/python"
```

### Rule: Historical variable names preserved but flagged

`EPGAT_PYTHON` (used in `frozen_protocol_benchmark.smk`, `label_scarcity_benchmark.smk`, `Figure4_graph_robustness.smk`) is a historical artifact. It actually points to the general runtime Python, not an EPGAT-specific interpreter. Do not rename without a standalone migration audit, but flag with comment:

```python
EPGAT_PYTHON = config["runtime"]["epgat_python_bin"]  # Historical name — points to general runtime Python
```

---

## 3. Wildcard Constraints

### Rule: Every wildcard MUST have a `wildcard_constraints` entry

Files using `{wildcards.xxx}` but lacking constraints:

| File | Wildcards used | Has constraints? |
|------|---------------|-----------------|
| `frozen_protocol_benchmark.smk` | protocol, model, feature, seed, data_key | ✅ Yes |
| `label_scarcity_benchmark.smk` | fraction_tag, model_variant, seed | ✅ Yes |
| `Figure4_graph_robustness.smk` | threshold, seed, graph_variant | ✅ Yes |
| `Figure2a/2b` | model, seed, feature | ✅ Yes |
| `Figure3a/3b/3c*` | model, seed, dim, variant | ✅ Yes |
| `Figure5_representation_mechanism.smk` | None (delegates to sub-workflows) | ✅ N/A |
| `Figure5d_feature_group_attribution.smk` | None | ✅ N/A |
| `fgraminearum_label_materialization.smk` | None | ✅ N/A |

### Pattern

```python
wildcard_constraints:
    protocol="|".join(PROTOCOLS),        # enum constraint
    model="|".join(MODEL_NAMES),         # enum constraint
    seed=r"\d+",                         # regex constraint
    threshold=r"\d+",                    # regex constraint
```

---

## 4. File Formatting (snakefmt)

### Tool: `snakefmt` (v2.0.3 tested)

Run before commit:
```bash
snakefmt --no-format-shell workflow/*.smk  # if shfmt unavailable
snakefmt workflow/*.smk                     # with shfmt installed
```

### Current status (2026-07-12)

- 15 of 19 `.smk` files would be reformatted (mostly line-wrapping + blank-line additions)
- 4 files already conform
- `shfmt` binary not available — shell directives formatted by hand for now

### Key rules snakefmt enforces

1. Blank line after `configfile:` directive
2. Line length: long lines broken with parentheses (`88` char soft limit)
3. Indentation: 4 spaces
4. List/dict trailing commas for multi-line
5. Two blank lines between top-level rules

### Install shfmt

```bash
conda install -n bioinfo -c conda-forge shfmt
# or
go install mvdan.cc/sh/v3/cmd/shfmt@latest
```

---

## 5. Workflow Structure Patterns

### Self-contained workflows

Each `workflow/Figure*.smk` must be independently executable:

```python
configfile: "configs/frozen_protocol.yaml"

# ... variable setup ...

rule freeze_protocol:  # Each workflow defines its own freeze rule
    ...

rule all:
    default_target: True
    input: ...
```

**Implication**: `freeze_protocol` is duplicated across 10+ workflows. This is intentional for independence but means a combined Snakefile would have rule conflicts. Do NOT create a monolithic Snakefile without resolving freeze_protocol deduplication first.

### Shared helpers via `include:`

```python
# workflow/frozen_protocol_benchmark.smk
include: "plots.smk"
```

```python
# workflow/Figure5a_hidden_umap_error_transition.smk (48-byte stub)
include: "../workflow/Figure5_representation_mechanism.smk"
```

**Note**: The 48-byte Figure5 stub files are valid proxy entry points. Document their purpose clearly.

### Output directory isolation

Every experiment must write to a unique output root:

```python
OUTPUT_ROOT = config["paths"].get("figure1_benchmark_output_root", "outputs/Figure1")
RESULTS_ROOT = config["paths"].get("figure1_results_root", "results/Figure1")
```

Never mix outputs from different experiments in the same directory.

---

## 6. Shell Directive Patterns

### Environment setup (standard block)

```bash
mkdir -p "{params.output_dir}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}"
export MPLBACKEND=Agg
export MPLCONFIGDIR="{MPLCONFIGDIR}"
export XDG_CACHE_HOME="{XDG_CACHE_HOME}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONPATH="${{PYTHONPATH:-.}}:."
```

This block is duplicated across ~13 rules. Consider extracting to a Snakemake `envvars:` directive or a shared shell snippet if Snakemake version supports it.

### Python module invocations

```bash
"{PYTHON_BIN}" -m src.module.entry_point \
  --config "{config_path}" \
  --protocol "{wildcards.protocol}" \
  > "{log}" 2>&1
```

Always use `{PYTHON_BIN}` (from config), not bare `python` or hardcoded paths.

---

## 7. Checklist for New Workflow Files

- [ ] Rule names follow `verb_object` lowercase convention
- [ ] All paths sourced from `config` dict, no hardcoded paths
- [ ] `wildcard_constraints` entry for every wildcard used
- [ ] `configfile:` directive as first non-comment line
- [ ] Blank line after `configfile:` line
- [ ] `rule all:` with `default_target: True`
- [ ] Each rule has `log:` for stdout/stderr capture
- [ ] `mkdir -p` before writing to output dirs
- [ ] Python invocation uses `{PYTHON_BIN}` from config
- [ ] Output root uses versioned/named subdirectory
- [ ] `snakefmt --check` passes before commit
- [ ] `snakemake -n --quiet all` succeeds as dry-run
