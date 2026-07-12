configfile: "configs/fgraminearum_label_materialization.yaml"

CONFIG_PATH = "configs/fgraminearum_label_materialization.yaml"
PYTHON_BIN = "/home/jiehuang/anaconda3/bin/python"

OLDLABEL_DIR = config["paths"]["oldlabel_dir"]
NEWLABEL_DIR = config["paths"]["newlabel_dir"]
BRIDGE_DIR = config["paths"]["bridge_dir"]
PROTOCOL_SOURCE_DIR = config["paths"]["protocol_source_dir"]
COMPARISON_TSV = config["paths"]["comparison_tsv"]
COMPARISON_MD = config["paths"]["comparison_md"]

BRIDGE_OUTPUTS = [
    BRIDGE_DIR + "/protein_to_canonical_bridge.tsv",
    BRIDGE_DIR + "/source_to_canonical_mapping.tsv",
    BRIDGE_DIR + "/high_confidence_yeast_transfer_candidates.tsv",
    BRIDGE_DIR + "/unresolved_high_confidence_ids.tsv",
    BRIDGE_DIR + "/bridge_summary.tsv",
    BRIDGE_DIR + "/bridge_summary.md",
    BRIDGE_DIR + "/bridge_source_manifest.tsv",
]

ALL_OUTPUTS = [
    OLDLABEL_DIR + "/labels.tsv",
    OLDLABEL_DIR + "/positive_genes.tsv",
    OLDLABEL_DIR + "/negative_genes.tsv",
    OLDLABEL_DIR + "/split.tsv",
    OLDLABEL_DIR + "/summary.tsv",
    OLDLABEL_DIR + "/summary.md",
    OLDLABEL_DIR + "/source_manifest.tsv",
    OLDLABEL_DIR + "/label_construction_audit.tsv",
    OLDLABEL_DIR + "/build_metadata.json",
    NEWLABEL_DIR + "/labels.tsv",
    NEWLABEL_DIR + "/positive_genes.tsv",
    NEWLABEL_DIR + "/negative_genes.tsv",
    NEWLABEL_DIR + "/split.tsv",
    NEWLABEL_DIR + "/summary.tsv",
    NEWLABEL_DIR + "/summary.md",
    NEWLABEL_DIR + "/source_manifest.tsv",
    NEWLABEL_DIR + "/label_construction_audit.tsv",
    NEWLABEL_DIR + "/build_metadata.json",
    COMPARISON_TSV,
    COMPARISON_MD,
]


rule build_fgraminearum_newlabel_bridge:
    input:
        config["paths"]["proteome_manifest"],
        config["paths"]["anchor_ncbi_protein_fasta"],
        config["paths"]["ph1_legacy_protein_fasta"],
        config["paths"]["ph1_legacy_mapping_tab"],
        config["paths"]["ph1_unified_id_map"],
        config["paths"]["yeast_transfer_table"],
        config["paths"]["master_evidence_mirror"],
        "data/processed/PPI/fgraminearum/string_id_mapping.tsv",
        "data/processed/EXP/fgraminearum/exp_id_mapping.tsv",
        "data/processed/LC/fgraminearum/subloc_id_mapping.tsv",
    output:
        BRIDGE_OUTPUTS
    log:
        BRIDGE_DIR + "/build_bridge.log",
    shell:
        r"""
        mkdir -p "{BRIDGE_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.data.build_fgraminearum_newlabel_bridge \
          --config "{CONFIG_PATH}" \
          > "{log}" 2>&1
        """


rule prepare_fgraminearum_label_materialization_sources:
    input:
        BRIDGE_OUTPUTS,
        config["paths"]["old_gene_list"],
        config["paths"]["yeast_transfer_table"],
        config["paths"]["master_evidence_mirror"],
    output:
        config["paths"]["old440_mapping_audit"],
        config["paths"]["old440_label_summary"],
        config["paths"]["lethal_positive_gene_list"],
    log:
        PROTOCOL_SOURCE_DIR + "/prepare_sources.log",
    shell:
        r"""
        mkdir -p "{PROTOCOL_SOURCE_DIR}"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.data.prepare_fgraminearum_label_materialization_sources \
          --config "{CONFIG_PATH}" \
          > "{log}" 2>&1
        """


rule materialize_fgraminearum_label_regimes:
    input:
        BRIDGE_OUTPUTS,
        config["paths"]["old_gene_list"],
        config["paths"]["old440_mapping_audit"],
        config["paths"]["old440_label_summary"],
        config["paths"]["lethal_positive_gene_list"],
        config["paths"]["yeast_transfer_table"],
        config["paths"]["master_evidence_mirror"],
    output:
        ALL_OUTPUTS
    log:
        "data/processed/essential_gene/fgraminearum/materialization.log"
    shell:
        r"""
        mkdir -p "data/processed/essential_gene/fgraminearum"
        export PYTHONPATH="${{PYTHONPATH:-.}}:."
        "{PYTHON_BIN}" -m src.data.materialize_fgraminearum_label_regimes \
          --config "{CONFIG_PATH}" \
          > "{log}" 2>&1
        """


rule all:
    default_target: True
    input:
        ALL_OUTPUTS
