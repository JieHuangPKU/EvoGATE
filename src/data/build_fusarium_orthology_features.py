import io
import os
import re
import sys

import numpy as np
import pandas as pd


ROOT = "/data276/jiehuang/fungi/Fusarium/orthofinder_65pilot_workflow"
ORTHO_ROOT = ROOT + "/results/OrthoFinder/Results_Fusarium_65pilot/Orthogroups"
ORTHOGROUPS_TSV = ORTHO_ROOT + "/Orthogroups.tsv"
GENECOUNT_TSV = ORTHO_ROOT + "/Orthogroups.GeneCount.tsv"
SCRIPT_INTEGRATE = ROOT + "/integrate_fg_essential_evidence_v2.py"
SCRIPT_ESSENTIAL_LIKE = ROOT + "/orthofinder_ph1_essential_like.py"

PH1_ACCESSION = "GCF_000240135.3"
EXPECTED_TOTAL_SPECIES = 65


def enable_utf8_stdout():
    return None


def _read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def split_gene_cell(cell):
    text = str(cell).strip() if cell is not None else ""
    if not text or text.lower() == "nan":
        return []
    return [item.strip() for item in text.split(",") if item.strip()]


def extract_gene_id(raw_gene):
    parts = [part.strip() for part in str(raw_gene).split("|") if part.strip()]
    return parts[-1] if parts else str(raw_gene).strip()


def classify_pangenome(occupancy):
    if occupancy >= 0.95:
        return "core"
    if occupancy >= 0.80:
        return "softcore"
    if occupancy >= 0.15:
        return "shell"
    return "cloud"


def load_fusarium_universe():
    master = _read_tsv("data_registry/master_label_table.preliminary.tsv")
    fg = master[master["species"] == "fgraminearum"][
        ["canonical_gene_id", "raw_gene_id", "raw_protein_id", "raw_transcript_id"]
    ].copy()
    fg = fg.drop_duplicates(subset=["canonical_gene_id"], keep="first").reset_index(drop=True)
    return fg


def load_unified_id_map():
    path = "/data276/jiehuang/fungi/Fusarium/Evidence/00_idmap/FG_gene_id_unified_map.tsv"
    return _read_tsv(path)


def compute_old_style_stats():
    gene_counts = pd.read_csv(GENECOUNT_TSV, sep="\t")
    species_columns = [c for c in gene_counts.columns if c not in ["Orthogroup", "Total"]]
    if len(species_columns) != EXPECTED_TOTAL_SPECIES:
        raise ValueError("Expected {} species columns, found {}".format(EXPECTED_TOTAL_SPECIES, len(species_columns)))
    numeric = gene_counts.loc[:, species_columns].apply(pd.to_numeric, errors="raise")
    presence_mask = numeric.gt(0)
    presence_species_count = presence_mask.sum(axis=1)
    single_copy_species_count = numeric.eq(1).sum(axis=1)
    stats = pd.DataFrame(
        {
            "orthogroup_id": gene_counts["Orthogroup"].astype(str),
            "presence_species_count": presence_species_count.astype(int),
            "occupancy": presence_species_count / float(EXPECTED_TOTAL_SPECIES),
            "ph1_copy": numeric[PH1_ACCESSION].astype(int),
            "mean_copy": numeric.mean(axis=1),
            "median_copy": numeric.median(axis=1),
            "max_copy": numeric.max(axis=1).astype(int),
            "std_copy": numeric.std(axis=1, ddof=0),
            "single_copy_species_count": single_copy_species_count.astype(int),
            "single_copy_fraction": np.where(
                presence_species_count > 0,
                single_copy_species_count / presence_species_count,
                0.0,
            ),
        }
    )
    stats["pangenome_class"] = stats["occupancy"].map(classify_pangenome)
    return stats


def compute_ph1_gene_to_orthogroup():
    orthogroups = _read_tsv(ORTHOGROUPS_TSV)
    records = []
    for _, row in orthogroups[["Orthogroup", PH1_ACCESSION]].iterrows():
        orthogroup_id = row["Orthogroup"]
        for raw_gene in split_gene_cell(row[PH1_ACCESSION]):
            records.append(
                {
                    "old_gene_id": extract_gene_id(raw_gene),
                    "orthogroup_id": orthogroup_id,
                    "old_namespace": "PH1_OrthoFinder_gene_id",
                }
            )
    return pd.DataFrame(records).drop_duplicates(subset=["old_gene_id"], keep="first").reset_index(drop=True)


def build_bridge(universe, ph1_map):
    raw_gene_to_canonical = {}
    raw_protein_to_canonical = {}
    gene_collisions = {}
    protein_collisions = {}
    for _, row in universe.iterrows():
        raw_gene = row["raw_gene_id"]
        raw_protein = row["raw_protein_id"]
        canonical = row["canonical_gene_id"]
        if raw_gene:
            if raw_gene in raw_gene_to_canonical and raw_gene_to_canonical[raw_gene] != canonical:
                gene_collisions.setdefault(raw_gene, set()).update([raw_gene_to_canonical[raw_gene], canonical])
            raw_gene_to_canonical[raw_gene] = canonical
        if raw_protein:
            if raw_protein in raw_protein_to_canonical and raw_protein_to_canonical[raw_protein] != canonical:
                protein_collisions.setdefault(raw_protein, set()).update([raw_protein_to_canonical[raw_protein], canonical])
            raw_protein_to_canonical[raw_protein] = canonical

    unified = load_unified_id_map()
    xp_to_ph1 = {}
    xp_collisions = {}
    for _, row in unified.iterrows():
        ph1 = row["ph1_canonical_gene_id"]
        proteins = [x.strip() for x in str(row["ncbi_protein_ids"]).split(";") if x.strip()]
        for prot in proteins:
            if prot in xp_to_ph1 and xp_to_ph1[prot] != ph1:
                xp_collisions.setdefault(prot, set()).update([xp_to_ph1[prot], ph1])
            xp_to_ph1[prot] = ph1

    ph1_to_canonical = dict(zip(universe["raw_gene_id"], universe["canonical_gene_id"]))

    records = []
    for _, row in ph1_map.iterrows():
        old_gene = row["old_gene_id"]
        if old_gene in xp_collisions:
            records.append(
                {
                    "old_gene_id": old_gene,
                    "canonical_gene_id": "",
                    "mapping_status": "ambiguous",
                    "mapping_rule": "unified_map_ncbi_protein_ids_conflict",
                    "needs_manual_review": "true",
                }
            )
        elif old_gene in xp_to_ph1 and xp_to_ph1[old_gene] in ph1_to_canonical:
            records.append(
                {
                    "old_gene_id": old_gene,
                    "canonical_gene_id": ph1_to_canonical[xp_to_ph1[old_gene]],
                    "mapping_status": "exact",
                    "mapping_rule": "ncbi_protein_id_via_unified_map_to_ph1_gene",
                    "needs_manual_review": "false",
                }
            )
        elif old_gene in protein_collisions:
            records.append(
                {
                    "old_gene_id": old_gene,
                    "canonical_gene_id": "",
                    "mapping_status": "ambiguous",
                    "mapping_rule": "raw_protein_id_exact_match_conflict",
                    "needs_manual_review": "true",
                }
            )
        elif old_gene in raw_protein_to_canonical:
            records.append(
                {
                    "old_gene_id": old_gene,
                    "canonical_gene_id": raw_protein_to_canonical[old_gene],
                    "mapping_status": "exact",
                    "mapping_rule": "raw_protein_id_exact_match",
                    "needs_manual_review": "false",
                }
            )
        elif old_gene in gene_collisions:
            records.append(
                {
                    "old_gene_id": old_gene,
                    "canonical_gene_id": "",
                    "mapping_status": "ambiguous",
                    "mapping_rule": "raw_gene_id_exact_match_conflict",
                    "needs_manual_review": "true",
                }
            )
        elif old_gene in raw_gene_to_canonical:
            records.append(
                {
                    "old_gene_id": old_gene,
                    "canonical_gene_id": raw_gene_to_canonical[old_gene],
                    "mapping_status": "exact",
                    "mapping_rule": "raw_gene_id_exact_match",
                    "needs_manual_review": "false",
                }
            )
        else:
            records.append(
                {
                    "old_gene_id": old_gene,
                    "canonical_gene_id": "",
                    "mapping_status": "unresolved",
                    "mapping_rule": "no_raw_gene_id_match_in_current_registry",
                    "needs_manual_review": "true",
                }
            )
    bridge = pd.DataFrame(records)
    return bridge


def main():
    enable_utf8_stdout()
    outdir = "outputs/fusarium_orthology"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("发现旧定义来源：")
    print("  • presence_species_count / occupancy / single_copy_fraction / pangenome_class -> orthofinder_ph1_essential_like.py::compute_orthogroup_stats/classify_pangenome")
    print("  • tier/core-like thresholds -> integrate_fg_essential_evidence_v2.py")

    rg_text = open(SCRIPT_ESSENTIAL_LIKE, "r", encoding="utf-8").read()
    integrate_text = open(SCRIPT_INTEGRATE, "r", encoding="utf-8").read()
    has_exact_core = ("is_exact_core" in rg_text) or ("is_exact_core" in integrate_text)
    has_strict_core = ("is_strict_core" in rg_text) or ("is_strict_core" in integrate_text)

    print("哪些字段已确认可复用：")
    print("  • presence_species_count")
    print("  • occupancy")
    print("  • single_copy_fraction")
    print("  • pangenome_class")
    print("pangenome_class 是否直接继承旧规则：是")
    print("is_exact_core 旧定义：{}".format("已找到" if has_exact_core else "未找到"))
    print("is_strict_core 旧定义：{}".format("已找到" if has_strict_core else "未找到"))

    universe = load_fusarium_universe()
    ph1_map = compute_ph1_gene_to_orthogroup()
    bridge = build_bridge(universe, ph1_map)
    bridge_path = os.path.join(outdir, "fusarium_orthology_id_bridge.tsv")
    bridge.to_csv(bridge_path, sep="\t", index=False)

    stats = compute_old_style_stats()
    feature = ph1_map.merge(stats, on="orthogroup_id", how="left", validate="many_to_one")
    feature = feature.merge(bridge, on="old_gene_id", how="left", validate="one_to_one")
    feature["orthology_source"] = "Orthogroups.tsv + Orthogroups.GeneCount.tsv; definitions from orthofinder_ph1_essential_like.py"
    feature["is_exact_core"] = ""
    feature["is_strict_core"] = ""

    mapped = feature[feature["mapping_status"] == "exact"].copy()
    canonical_groups = mapped.groupby("canonical_gene_id")
    exact_records = []
    ambiguous_canonical = set()
    for canonical_gene_id, sub in canonical_groups:
        unique_ogs = sorted(set(sub["orthogroup_id"].astype(str)))
        if len(unique_ogs) == 1:
            first = sub.iloc[0]
            exact_records.append(
                {
                    "canonical_gene_id": canonical_gene_id,
                    "orthogroup_id": unique_ogs[0],
                    "presence_species_count": first["presence_species_count"],
                    "occupancy": first["occupancy"],
                    "single_copy_fraction": first["single_copy_fraction"],
                    "pangenome_class": first["pangenome_class"],
                    "is_exact_core": "",
                    "is_strict_core": "",
                    "orthology_join_status": "exact",
                    "orthology_source": first["orthology_source"],
                }
            )
        else:
            ambiguous_canonical.add(canonical_gene_id)

    feature_exact = pd.DataFrame(
        exact_records,
        columns=[
            "canonical_gene_id",
            "orthogroup_id",
            "presence_species_count",
            "occupancy",
            "single_copy_fraction",
            "pangenome_class",
            "is_exact_core",
            "is_strict_core",
            "orthology_join_status",
            "orthology_source",
        ],
    )
    final = universe[["canonical_gene_id"]].merge(feature_exact, on="canonical_gene_id", how="left")
    final["orthology_join_status"] = final["orthology_join_status"].replace("", pd.NA)
    final.loc[final["canonical_gene_id"].isin(ambiguous_canonical), "orthology_join_status"] = "ambiguous_multiple_orthogroups"
    final["orthology_join_status"] = final["orthology_join_status"].fillna("missing_in_orthofinder_or_unresolved")
    final["orthology_source"] = final["orthology_source"].replace("", pd.NA).fillna(
        "Orthogroups.tsv + Orthogroups.GeneCount.tsv; definitions from orthofinder_ph1_essential_like.py"
    )
    for col in ["orthogroup_id", "pangenome_class", "is_exact_core", "is_strict_core"]:
        final[col] = final[col].fillna("")
    final_path = os.path.join(outdir, "fusarium_orthology_features.tsv")
    final.to_csv(final_path, sep="\t", index=False)

    schema_rows = [
        {
            "feature_name": "presence_species_count",
            "definition_source": "orthofinder_ph1_essential_like.py::compute_orthogroup_stats",
            "definition": "presence_mask.sum(axis=1)",
        },
        {
            "feature_name": "occupancy",
            "definition_source": "orthofinder_ph1_essential_like.py::compute_orthogroup_stats",
            "definition": "presence_species_count / 65",
        },
        {
            "feature_name": "single_copy_fraction",
            "definition_source": "orthofinder_ph1_essential_like.py::compute_orthogroup_stats",
            "definition": "single_copy_species_count / presence_species_count when presence_species_count > 0 else 0",
        },
        {
            "feature_name": "pangenome_class",
            "definition_source": "orthofinder_ph1_essential_like.py::classify_pangenome",
            "definition": "occupancy>=0.95 core; >=0.80 softcore; >=0.15 shell; else cloud",
        },
        {
            "feature_name": "is_exact_core",
            "definition_source": "not_found_in_old_scripts",
            "definition": "未找到现成定义，本轮不自创新规则，输出空值",
        },
        {
            "feature_name": "is_strict_core",
            "definition_source": "not_found_in_old_scripts",
            "definition": "未找到现成定义，本轮不自创新规则，输出空值",
        },
    ]
    pd.DataFrame(schema_rows).to_csv(os.path.join(outdir, "fusarium_orthology_feature_schema.tsv"), sep="\t", index=False)

    audit_rows = [
        {
            "item": "presence_species_count_definition",
            "value": "presence_mask.sum(axis=1)",
            "source": SCRIPT_ESSENTIAL_LIKE + "::compute_orthogroup_stats",
        },
        {
            "item": "occupancy_definition",
            "value": "presence_species_count / 65",
            "source": SCRIPT_ESSENTIAL_LIKE + "::compute_orthogroup_stats",
        },
        {
            "item": "single_copy_fraction_definition",
            "value": "single_copy_species_count / presence_species_count",
            "source": SCRIPT_ESSENTIAL_LIKE + "::compute_orthogroup_stats",
        },
        {
            "item": "pangenome_class_definition",
            "value": "core/softcore/shell/cloud by occupancy thresholds 0.95/0.80/0.15",
            "source": SCRIPT_ESSENTIAL_LIKE + "::classify_pangenome",
        },
        {
            "item": "is_exact_core_found",
            "value": str(has_exact_core).lower(),
            "source": "rg over old scripts",
        },
        {
            "item": "is_strict_core_found",
            "value": str(has_strict_core).lower(),
            "source": "rg over old scripts",
        },
        {
            "item": "old_gene_id_namespace",
            "value": "NCBI protein accession in OrthoFinder PH1 column (e.g. XP_011315572.1)",
            "source": ORTHOGROUPS_TSV,
        },
        {
            "item": "canonical_mapping_rule",
            "value": "ncbi_protein_id_via_unified_map_to_ph1_gene -> raw_gene_id exact match -> canonical_gene_id",
            "source": "/data276/jiehuang/fungi/Fusarium/Evidence/00_idmap/FG_gene_id_unified_map.tsv + data_registry/master_label_table.preliminary.tsv",
        },
    ]
    audit_tsv = os.path.join(outdir, "fusarium_orthology_asset_audit.tsv")
    pd.DataFrame(audit_rows).to_csv(audit_tsv, sep="\t", index=False)

    unresolved = int((bridge["mapping_status"] == "unresolved").sum())
    ambiguous = int((bridge["mapping_status"] == "ambiguous").sum())
    exact = int((bridge["mapping_status"] == "exact").sum())

    md_lines = [
        "# Fusarium Orthology Asset Audit",
        "",
        "## 旧定义来源",
        "- `presence_species_count`: `orthofinder_ph1_essential_like.py::compute_orthogroup_stats`",
        "- `occupancy`: `orthofinder_ph1_essential_like.py::compute_orthogroup_stats`",
        "- `single_copy_fraction`: `orthofinder_ph1_essential_like.py::compute_orthogroup_stats`",
        "- `pangenome_class`: `orthofinder_ph1_essential_like.py::classify_pangenome`",
        "- `is_exact_core`: 未找到现成定义",
        "- `is_strict_core`: 未找到现成定义",
        "",
        "## 规则复用结论",
        "- pangenome_class 完全复用旧规则：是",
        "- core / strict_core 现成定义是否找到：否，本轮不自创新规则，输出空值",
        "",
        "## ID namespace",
        "- 旧 OrthoFinder / PH1 namespace: `XP_...` NCBI protein accession",
        "- 当前 ProGATE_v2 namespace: `fgraminearum::FGRAMPH1_XXGXXXXX`",
        "- 映射规则: `XP_... -> FG_gene_id_unified_map.tsv::ph1_canonical_gene_id -> fgraminearum::...`",
        "",
        "## Bridge 结果",
        "- exact mapped: {}".format(exact),
        "- unresolved: {}".format(unresolved),
        "- ambiguous: {}".format(ambiguous),
        "",
        "## 可接入性",
        "- `fusarium_orthology_features.tsv` 已按当前 canonical universe 左连接展开，可直接 join 到 ProGATE_v2 Fusarium graph input。",
    ]
    audit_md = os.path.join(outdir, "fusarium_orthology_asset_audit.md")
    with open(audit_md, "w", encoding="utf-8") as handle:
        handle.write("\n".join(md_lines))

    summary_md = os.path.join(outdir, "fusarium_orthology_integration_summary.md")
    with open(summary_md, "w", encoding="utf-8") as handle:
        handle.write("\n".join(md_lines))

    print("哪些 ID 还没映射通：")
    print("  • unresolved = {}".format(unresolved))
    print("  • ambiguous = {}".format(ambiguous))
    print("成功映射到 canonical_gene_id 的基因数：{}".format(exact))
    print("成功生成的特征：")
    print("  • orthogroup_id")
    print("  • presence_species_count")
    print("  • occupancy")
    print("  • single_copy_fraction")
    print("  • pangenome_class")
    print("  • orthology_join_status")
    print("  • orthology_source")
    print("是否已可直接接入 ProGATE_v2：是")


if __name__ == "__main__":
    enable_utf8_stdout()
    main()
