from __future__ import annotations

from dataclasses import dataclass, asdict


@dataclass
class SupportGraphRegistryRecord:
    species: str
    graph_completeness_class: str
    node_table_path: str
    edge_table_path: str
    node_feature_manifest_path: str
    label_manifest_path: str
    edge_weight_default: float
    use_in_graph_training: bool
    notes: str

    def to_dict(self) -> dict:
        return asdict(self)
