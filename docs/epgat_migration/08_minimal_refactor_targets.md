# Minimal Refactor Targets

## Refactor First

### `utils/utils.py`

Reason:
- defines the real legacy dataset contract
- currently hides core assumptions about labels, features, and PPI filtering

Action:
- extract into explicit dataset builder plus schema manifest

### `runners/tools.py`

Reason:
- central bridge for CLI flags, feature assembly, embedding injection, and graph preprocessing

Action:
- split into:
  - argument compatibility layer
  - feature assembly layer
  - graph preprocessing layer

### `runners/run_gat.py`

Reason:
- closest to original EPGAT core but now mixed with expanded evaluation side effects

Action:
- preserve model behavior
- move output/export logic out

### Fusarium Legacy Feature Inputs

Targets:
- `data/essential_genes/fgraminearum/EssentialGenes/ogee.csv`
- `data/essential_genes/fgraminearum/PPI/STRING/string.csv`
- `data/essential_genes/fgraminearum/Expression/profile.csv`
- `data/essential_genes/fgraminearum/Orthologs/orthologs.csv`
- `data/essential_genes/fgraminearum/SubLocalizations/subloc.csv`

Reason:
- these need canonical adapters before any deeper model work

## Preserve For Now

### `models/gat/gat_pytorch.py`

Reason:
- active core implementation
- preserve behavior during first migration stage

### `models/gcn/gcn_pytorch.py`

Reason:
- active implementation and cleaner than some surrounding runner logic

### `models/gin/gin_pytorch.py`

Reason:
- active implementation

### `models/graphsage/graphsage_adapter.py`

Reason:
- active adapter logic likely still needed for legacy replay

## Do Not Refactor Early

### `utils/prepare_data.py`

Reason:
- not suitable as trusted runtime logic
- contains hard-coded local paths and post-hoc metric manipulation

Recommendation:
- keep as historical reference only

### `runners/backup/` And `runners/20250714/`

Reason:
- stale copies and dated snapshots

Recommendation:
- exclude from migration implementation

### `results/Third/`

Reason:
- duplicate cohort outputs rather than code logic

Recommendation:
- keep only as provenance reference
