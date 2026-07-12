import os

import pandas as pd


def build_support_prior_matrices(best_model_name, oof_scores_path, output_dir, species_list, node_table_template):
    scores = pd.read_csv(oof_scores_path, sep="\t", dtype=str).fillna("")
    best = scores[scores["model_name"] == best_model_name].copy()
    written = []
    for species in species_list:
        nodes = pd.read_csv(node_table_template.format(species=species), sep="\t", dtype=str).fillna("")
        sub = best[best["species"] == species][["canonical_gene_id", "prior_score"]].copy()
        sub["prior_score"] = pd.to_numeric(sub["prior_score"], errors="coerce").fillna(0.0)
        matrix = nodes[["species", "canonical_gene_id"]].merge(sub, on="canonical_gene_id", how="left")
        matrix["prior_score"] = pd.to_numeric(matrix["prior_score"], errors="coerce").fillna(0.0)
        matrix["has_prior_score"] = (matrix["prior_score"] != 0).astype(float)
        matrix["prior_missing_mask"] = (matrix["has_prior_score"] <= 0).astype(float)
        path = os.path.join(output_dir, "support_prior_matrix_{}.tsv".format(species))
        matrix.to_csv(path, sep="\t", index=False)
        written.append(path)
    return written


def write_fusarium_prior(best_model_name, fusarium_scores_path, output_path):
    scores = pd.read_csv(fusarium_scores_path, sep="\t", dtype=str).fillna("")
    best = scores[scores["model_name"] == best_model_name].copy()
    best.to_csv(output_path, sep="\t", index=False)
    return output_path
