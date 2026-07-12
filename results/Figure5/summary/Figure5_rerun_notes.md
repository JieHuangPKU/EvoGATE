# Figure5 Rerun Notes

1. Run `scripts/run_Figure5_fusarium_graphsage_bio_esm2_representation_mechanism.sh -j 48` to rebuild the full Figure5 module.
2. Run `scripts/run_Figure5a_fusarium_graphsage_bio_esm2_hidden_umap_error_transition.sh -j 48`, `scripts/run_Figure5b_fusarium_graphsage_bio_esm2_hidden_quant_summary.sh -j 48`, or `scripts/run_Figure5c_fusarium_graphsage_bio_esm2_input_vs_hidden_compare.sh -j 48` for panel-specific rebuilds.
3. The legacy export-style Figure5 representation plots are no longer part of the active Figure5 DAG.
