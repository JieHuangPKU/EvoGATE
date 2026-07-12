"""
Original-compatible GraphSAGE for the archived replay benchmark.

This implementation keeps the old runner contract while fixing the
aggregator_type bug. We implement the mean/pool aggregators directly
so the behavior is stable across the old PyG environment.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd

try:
    from torch_scatter import scatter_max, scatter_mean
except ImportError as exc:
    raise ImportError("GraphSAGE requires torch_scatter in the EPGAT environment") from exc

from src.models.multimodal_gated_fusion import MultimodalGatedFusion


class _OriginalCompatibleSAGELayer(nn.Module):
    def __init__(self, in_feats, out_feats, aggregator_type="mean"):
        super().__init__()
        normalized = str(aggregator_type).strip().lower()
        if normalized not in {"mean", "pool"}:
            raise ValueError(f"Unsupported GraphSAGE aggregator_type '{aggregator_type}'")
        self.aggregator_type = normalized
        self.lin = nn.Linear(int(in_feats) * 2, int(out_feats))
        self.pool_linear = nn.Linear(int(in_feats), int(in_feats)) if self.aggregator_type == "pool" else None

    def forward(self, x, edge_index):
        src = edge_index[0]
        dst = edge_index[1]
        num_nodes = x.shape[0]
        if self.aggregator_type == "mean":
            neighbor_agg = scatter_mean(x[src], dst, dim=0, dim_size=num_nodes)
        else:
            pooled = F.relu(self.pool_linear(x[src]))
            neighbor_agg, _ = scatter_max(pooled, dst, dim=0, dim_size=num_nodes)
            neighbor_agg[torch.isinf(neighbor_agg)] = 0.0
        combined = torch.cat([x, neighbor_agg], dim=-1)
        return self.lin(combined)


class EPGATOriginalSAGE(nn.Module):
    def __init__(self, in_feats, n_hidden=64, n_layers=2, dropout=0.5, aggregator_type="mean"):
        super().__init__()
        self.dropout = dropout
        self.aggregator_type = str(aggregator_type).strip().lower()
        self.layers = nn.ModuleList()
        hidden_in = int(in_feats)
        hidden_layers = max(int(n_layers) - 1, 1)
        for _ in range(hidden_layers):
            self.layers.append(_OriginalCompatibleSAGELayer(hidden_in, int(n_hidden), aggregator_type=self.aggregator_type))
            hidden_in = int(n_hidden)
        self.out_layer = _OriginalCompatibleSAGELayer(hidden_in, 1, aggregator_type=self.aggregator_type)

    def forward(self, x, edge_index, edge_weights=None):
        for layer in self.layers:
            x = layer(x, edge_index)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        x = self.out_layer(x, edge_index)
        return x


class EPGATOriginalSAGEWithFusion(nn.Module):
    def __init__(
        self,
        in_feats,
        fusion_partition,
        fusion_hidden_dim=256,
        fusion_dropout=0.2,
        fusion_mode="gated",
        n_hidden=64,
        n_layers=2,
        dropout=0.5,
        aggregator_type="mean",
    ):
        super().__init__()
        self.fusion_partition = dict(fusion_partition)
        self.fusion = MultimodalGatedFusion(
            omics_dim=int(self.fusion_partition["omics_dim"]),
            esm2_dim=int(self.fusion_partition["esm2_dim"]),
            hidden_dim=int(fusion_hidden_dim),
            dropout=float(fusion_dropout),
            fusion_mode=str(fusion_mode),
        )
        self.encoder = EPGATOriginalSAGE(
            in_feats=int(self.fusion.output_dim),
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout=dropout,
            aggregator_type=aggregator_type,
        )
        self._last_gate = None

    def _split_modalities(self, x):
        omics = x[:, int(self.fusion_partition["omics_start_col"]): int(self.fusion_partition["omics_end_col"]) + 1]
        esm2 = x[:, int(self.fusion_partition["esm2_start_col"]): int(self.fusion_partition["esm2_end_col"]) + 1]
        return omics, esm2

    def forward(self, x, edge_index, edge_weights=None):
        x_omics, x_esm2 = self._split_modalities(x)
        fused, gate = self.fusion(x_omics, x_esm2)
        self._last_gate = gate.detach()
        return self.encoder(fused, edge_index, edge_weights=edge_weights)

    def gate_statistics(self, x, node_manifest=None, split_manifest=None):
        with torch.no_grad():
            x_omics, x_esm2 = self._split_modalities(x)
            _, gate = self.fusion(x_omics, x_esm2)
        gate_cpu = gate.detach().cpu()
        gate_mean_by_node = gate_cpu.mean(dim=1).numpy()
        rows = [
            {"subset": "all_nodes", "metric": "gate_mean", "value": float(gate_cpu.mean().item())},
            {"subset": "all_nodes", "metric": "gate_median", "value": float(gate_cpu.median().item())},
            {"subset": "all_nodes", "metric": "gate_std", "value": float(gate_cpu.std(unbiased=False).item())},
        ]
        dim_mean = gate_cpu.mean(dim=0)
        dim_std = gate_cpu.std(dim=0, unbiased=False)
        rows.extend(
            {
                "subset": "all_nodes",
                "metric": "gate_dim_mean",
                "dimension": int(idx),
                "value": float(value.item()),
            }
            for idx, value in enumerate(dim_mean)
        )
        rows.extend(
            {
                "subset": "all_nodes",
                "metric": "gate_dim_std",
                "dimension": int(idx),
                "value": float(value.item()),
            }
            for idx, value in enumerate(dim_std)
        )
        if node_manifest is not None and split_manifest is not None and not split_manifest.empty:
            node_gate = node_manifest[["graph_gene_id"]].copy()
            node_gate["gate_mean"] = gate_mean_by_node
            manifest = split_manifest.copy().merge(node_gate, on="graph_gene_id", how="left")
            manifest["label"] = pd.to_numeric(manifest["label"], errors="coerce")
            for split_name in ["train", "val", "test"]:
                split_mask = manifest["split"].astype(str) == split_name
                if split_mask.any():
                    rows.append(
                        {
                            "subset": split_name,
                            "metric": "gate_mean",
                            "value": float(manifest.loc[split_mask, "gate_mean"].mean()),
                        }
                    )
            for label_value, label_name in [(1, "essential"), (0, "non_essential")]:
                label_mask = manifest["label"] == label_value
                if label_mask.any():
                    rows.append(
                        {
                            "subset": label_name,
                            "metric": "gate_mean",
                            "value": float(manifest.loc[label_mask, "gate_mean"].mean()),
                        }
                    )
        return pd.DataFrame(rows)
