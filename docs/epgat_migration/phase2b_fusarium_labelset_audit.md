# Phase 2B Fusarium Labelset Audit

Update: the file paths below are retained as historical comparison snapshots. They are not the
current mainline Fusarium label inputs.

## Label set A: old 440 set

- Historical archive positive source: `results/label_rebuild_experiments/old440/labels/positive_old440.tsv`.
- Historical archive negative source: `results/label_rebuild_experiments/old440/labels/negative_old440.tsv`.
- Final size: 439 positives and 10154 negatives.
- ID system: canonical gene IDs after explicit mapping from `gene_list.txt` source IDs.
- Canonical mapping rate: 1.0000.

## Label set B: new lethal + yeast-transfer set

- Historical archive positive source: `results/label_rebuild_experiments/labels/positive_set_P1.tsv`.
- Historical archive negative source: `results/label_rebuild_experiments/labels/negative_set.tsv`.
- Final size: 1096 positives and 10270 negatives.
- ID system: canonical gene IDs already resolved in the rebuilt label compare workflow.
- Canonical mapping rate: 1.0000 over the final emitted tables.

## ID compatibility

- Both label schemes can already be consumed in canonical gene space.
- The old 440 positives preserve source->canonical audit columns; the new set is already canonicalized.
- Current mainline materialized label assets are now maintained under `data/processed/essential_gene/fgraminearum/{oldlabel,newlabel}/`.

## Recommended comparison protocol

- Reuse the same `fgraminearum` benchmark graph/model path and swap only the label manifest.
- Run the comparison on the top 1-2 archived replay models first, or after the human Phase 2A old-result target is stabilized.
- Keep the fixed legacy feature contract unchanged so the label effect is isolated.
