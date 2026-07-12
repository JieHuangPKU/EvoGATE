# EvoGATE migration guide

_Non-destructive guidance for the future ProGATE_v2-to-EvoGATE migration._

---

## Scope

This guide records known migration work. It does not authorize path replacement, file movement, renaming, symlink creation, or result regeneration. The current documentation phase performs no migration.

## Identity mapping

| Item | Historical value | Current value | Policy |
|---|---|---|---|
| Project name | ProGATE / ProGATE_v2 | EvoGATE | Use EvoGATE in new documentation and claims |
| Method predecessor | EPGAT | EPGAT | Preserve name; never relabel as EvoGATE |
| External comparator/source | Bingo | Bingo | Preserve name and external status |
| Historical Linux root | `/home/jiehuang/software/fungi/ProGATE_v2` | `/DATA/software/bioinfo/fungi/EvoGATE` | Report; do not batch-replace |
| Historical macOS root | `/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE_v2` | Current repository root | Report; migrate producer-by-producer later |

## Known hard-coded path classes

### Shell wrappers

Several `scripts/run_*.sh` wrappers assign `PROJECT_ROOT=/home/jiehuang/software/fungi/ProGATE_v2`. Examples include Figure1, Figure2, Figure3, Figure4, label-scarcity, and label-materialization wrappers. These are **Historical / non-portable**.

### Runtime configuration

`configs/frozen_protocol.yaml` contains historical `mplconfigdir`, cache, Python environment, and local ESM2 model paths. These paths document the environment used but are not portable defaults.

### Data builders

Several builders under `src/data/` contain historical macOS or external Linux paths for raw sources, EPGAT assets, and Fusarium evidence. Existing processed artifacts may be usable while their exact raw rebuild remains blocked.

### Generated reports

Result summaries and migration documents contain historical absolute output paths. These should remain unchanged when they record provenance. New documentation should refer to repository-relative paths.

## Historical entry points

| Historical entry type | Current interpretation |
|---|---|
| ProGATE_v2 Shell wrapper with absolute root | Non-portable; not recommended |
| EPGAT replay runner | Historical comparison only |
| Old ranking-only Fusarium workflow | Historical; not mainline frozen protocol |
| Figure workflow with current repository-relative inputs | Current experiment entry, subject to missing outputs/environment |
| `workflow/frozen_protocol_benchmark.smk` | Current canonical benchmark workflow |
| `workflow/fgraminearum_label_materialization.smk` | Current canonical label workflow |

## Future migration principles

1. Recover or establish version control before structural migration.
2. Inventory every reference before changing a path.
3. Separate provenance text from executable configuration.
4. Replace one executable path class at a time and leave a validation record.
5. Keep historical artifacts immutable.
6. Do not redirect old paths with root-level symlinks.
7. Do not combine path migration with scientific protocol changes.
8. Preserve oldlabel and legacy replay as explicit, non-default branches.
9. Write new runs to versioned output and result directories.
10. Validate all consumers before deprecating an old entry.

## Proposed migration sequence

| Phase | Change | Required evidence |
|---|---|---|
| 1 | Establish Git baseline and authoritative file inventory | Repository history or signed baseline manifest |
| 2 | Define a tested environment | Locked Python/R/Snakemake dependencies |
| 3 | Create portable root-relative canonical commands | Dry-run and lightweight checks |
| 4 | Recover missing source and outputs | Checksums and provenance mapping |
| 5 | Migrate data builders from absolute raw paths | Per-modality rebuild comparison |
| 6 | Consolidate workflow entry points | Target and output equivalence audit |
| 7 | Archive deprecated wrappers | Confirmed no current consumers |

## Migration stop conditions

Stop and request review if a path change alters a sample count, label count, identifier mapping, split, seed, graph edge set, feature dimension, metric, or output destination. Such a change is scientific or behavioral, not merely a path migration.

## Current blockers

Migration validation is **Blocked** by empty Git metadata, missing environment locks, missing per-run `outputs/`, bytecode-only modules, and the missing yeast-transfer confidence producer. These blockers should be resolved in separate approved work.

