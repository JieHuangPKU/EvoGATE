# Fusarium Label Materialization Dependency Fix

## What Was Still Wrong

`workflow/fgraminearum_label_materialization.smk` still declared these historical files as direct inputs to `materialize_fgraminearum_label_regimes`:

- `results/label_rebuild_experiments/old440/labels/old440_mapping_audit.tsv`
- `results/label_rebuild_experiments/old440/labels/old440_label_summary.tsv`
- `results/label_rebuild_experiments/labels/lethal_positive_gene_list.tsv`

This was not only a stale Snakemake declaration problem. `src/data/materialize_fgraminearum_label_regimes.py` also read the same config-backed paths at runtime:

- `build_newlabel()` read `paths["lethal_positive_gene_list"]`
- `build_oldlabel()` read `paths["old440_mapping_audit"]`
- `build_oldlabel()` read `paths["old440_label_summary"]`

At that historical point, the workflow still had a real runtime dependency on `results/label_rebuild_experiments/...`.

## Why Dry-Run And Force-Run Behaved Differently

Before this fix, plain dry-run targeted `rule all`. Because the processed outputs under `data/processed/essential_gene/fgraminearum/` already existed, Snakemake reported:

- `Nothing to be done (all requested files are present and up to date).`

That did not mean the dependency graph was clean. It only meant Snakemake did not need to schedule `materialize_fgraminearum_label_regimes` for that invocation.

When `--forcerun materialize_fgraminearum_label_regimes` was added, Snakemake was required to include that rule in the DAG and validate its declared inputs. At that point the hidden inconsistency became visible:

- the rule still formally required the missing legacy files
- the Python entrypoint still read the same legacy files from config

So dry-run looked clean only because the target outputs already existed, while force-run exposed the unresolved upstream inputs.

## What Was Fixed

The workflow now protocolizes the missing source tables inside this workflow instead of reading them from `results/label_rebuild_experiments/...`.

Modified files:

- `workflow/fgraminearum_label_materialization.smk`
- `configs/fgraminearum_label_materialization.yaml`
- `src/data/prepare_fgraminearum_label_materialization_sources.py`

### New Protocolized Source Step

A new rule was added:

- `prepare_fgraminearum_label_materialization_sources`

It writes these protocolized inputs under:

- `data/interim/protocol_refactor/fgraminearum_label_materialization/old440_mapping_audit.tsv`
- `data/interim/protocol_refactor/fgraminearum_label_materialization/old440_label_summary.tsv`
- `data/interim/protocol_refactor/fgraminearum_label_materialization/lethal_positive_gene_list.tsv`

### Reconstructed Inputs

`old440_mapping_audit.tsv`

- rebuilt from the historical `gene_list.txt`
- mapped through repo-local canonical ID logic
- uses master evidence raw-ID mapping and the protocolized bridge outputs when needed

`old440_label_summary.tsv`

- rebuilt from the reconstructed old440 replay audit
- negative counts now reflect the protocolized replay against the current rebuilt negative pool
- this is no longer a borrowed summary table from the historical exploratory results directory

`lethal_positive_gene_list.tsv`

- rebuilt from `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv`
- selection rule is now explicit and reproducible:
  `phi-base_current.csv` + `evidence_term_raw == lethal` + `supports_gold_label == true`
- transcript-only raw IDs (`FGRAMPH1_01T...`) are excluded
- rows are then normalized into the final `fgraminearum::FGRAMPH1_*` ID space using repo-local source mappings and the protocolized bridge
- PHI lethal rows that still cannot be normalized into the final canonical space are excluded from the runtime provenance table rather than leaking `FGSG_*` IDs into processed labels

## Current Runtime Dependency Graph

The materialization workflow now reads these protocolized source files instead of the historical `results/...` files:

- `data/interim/protocol_refactor/fgraminearum_label_materialization/old440_mapping_audit.tsv`
- `data/interim/protocol_refactor/fgraminearum_label_materialization/old440_label_summary.tsv`
- `data/interim/protocol_refactor/fgraminearum_label_materialization/lethal_positive_gene_list.tsv`

It still depends on these non-results upstream artifacts:

- `/home/jiehuang/software/fungi/EPGAT/data/essential_genes/fgraminearum/EssentialGenes/gene_list.txt`
- `data/derived_labels/ph1_yeast_essential_ortholog_labels.tsv`
- `data/interim/protocol_refactor/master_evidence_table.preliminary.tsv`
- `data/processed/essential_gene/fgraminearum/bridge/*`

## Is `results/label_rebuild_experiments/...` Still In The Runtime Graph?

No.

It remains in the repository only as a historical archive / provenance directory and should not
be treated as a current mainline workflow input.

After the fix:

- Snakemake inputs for this workflow no longer reference `results/label_rebuild_experiments/...`
- the config used by this workflow no longer points to `results/label_rebuild_experiments/...`
- the Python modules used at runtime by this workflow no longer read `results/label_rebuild_experiments/...`
- the generated `source_manifest.tsv` files for both `oldlabel` and `newlabel` now point to protocolized `data/interim/...` sources instead

## Validation

Commands used:

```bash
XDG_CACHE_HOME=/home/jiehuang/software/fungi/ProGATE_v2/.cache \
MPLCONFIGDIR=/home/jiehuang/software/fungi/ProGATE_v2/.mplconfig \
snakemake -n --printshellcmds --verbose \
  -s workflow/fgraminearum_label_materialization.smk \
  --configfile configs/fgraminearum_label_materialization.yaml \
  --cores 1
```

```bash
./scripts/run_fgraminearum_label_materialization.sh 16 --forcerun materialize_fgraminearum_label_regimes
```

```bash
XDG_CACHE_HOME=/home/jiehuang/software/fungi/ProGATE_v2/.cache \
MPLCONFIGDIR=/home/jiehuang/software/fungi/ProGATE_v2/.mplconfig \
snakemake -n --printshellcmds --verbose \
  -s workflow/fgraminearum_label_materialization.smk \
  --configfile configs/fgraminearum_label_materialization.yaml \
  --cores 1
```

Observed behavior:

- pre-run dry-run: scheduled the new protocolized source-prep rule and then materialization
- real `--forcerun materialize_fgraminearum_label_regimes`: completed successfully
- post-run dry-run: `Nothing to be done (all requested files are present and up to date).`
