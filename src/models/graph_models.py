from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class GraphDependencyStatus:
    torch_available: bool
    torch_geometric_available: bool
    dgl_available: bool
    notes: list[str]


def inspect_graph_dependencies() -> GraphDependencyStatus:
    notes: list[str] = []
    torch_available = False
    torch_geometric_available = False
    dgl_available = False

    try:
        import torch  # noqa: F401

        torch_available = True
        notes.append("torch: available")
    except Exception as exc:  # pragma: no cover
        notes.append(f"torch: missing ({type(exc).__name__}: {exc})")

    try:
        import torch_geometric  # noqa: F401

        torch_geometric_available = True
        notes.append("torch_geometric: available")
    except Exception as exc:  # pragma: no cover
        notes.append(f"torch_geometric: missing ({type(exc).__name__}: {exc})")

    try:
        import dgl  # noqa: F401

        dgl_available = True
        notes.append("dgl: available")
    except Exception as exc:  # pragma: no cover
        notes.append(f"dgl: missing ({type(exc).__name__}: {exc})")

    return GraphDependencyStatus(
        torch_available=torch_available,
        torch_geometric_available=torch_geometric_available,
        dgl_available=dgl_available,
        notes=notes,
    )


def _require_torch() -> Any:
    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("PyTorch is required for graph models in ProGATE_v2.") from exc
    return torch, nn, F


class _GraphRuntimeMixin:
    backend_name = "torch"
    model_name = "graph_model"

    def encode(self, x, edge_index):
        return self.forward(x, edge_index)

    def decode_anchor_scores(self, embeddings, anchor_mask, temperature: float = 5.0):
        torch, _, F = _require_torch()
        if anchor_mask.sum().item() == 0:
            # Fall back to normalized embedding norm if there are no anchors.
            scores = embeddings.norm(dim=1)
            scores = (scores - scores.mean()) / scores.std().clamp(min=1e-6)
            return torch.sigmoid(scores)
        anchor_embeddings = embeddings[anchor_mask]
        anchor_centroid = F.normalize(anchor_embeddings.mean(dim=0, keepdim=True), dim=1)
        normalized_embeddings = F.normalize(embeddings, dim=1)
        cosine_scores = (normalized_embeddings @ anchor_centroid.T).squeeze(1)
        return torch.sigmoid(cosine_scores * temperature)

    def reconstruction_loss(self, embeddings, positive_edge_index):
        raise NotImplementedError

    def model_summary(self) -> dict[str, Any]:
        return {
            "backend_name": self.backend_name,
            "model_name": self.model_name,
            "class_name": self.__class__.__name__,
            "parameter_count": int(sum(parameter.numel() for parameter in self.parameters())),
        }


def _build_torch_graphsage_model(
    input_dim: int,
    hidden_dim: int = 64,
    output_dim: int = 1,
    num_layers: int = 2,
    dropout: float = 0.2,
):
    torch, nn, _ = _require_torch()

    class GraphSAGELayer(nn.Module):
        def __init__(self, in_dim: int, out_dim: int) -> None:
            super().__init__()
            self.self_linear = nn.Linear(in_dim, out_dim)
            self.neighbor_linear = nn.Linear(in_dim, out_dim)

        def forward(self, x, edge_index):
            num_nodes = x.shape[0]
            src, dst = edge_index[0], edge_index[1]
            aggregated = torch.zeros_like(x)
            aggregated.index_add_(0, dst, x[src])
            degree = torch.zeros(num_nodes, device=x.device, dtype=x.dtype)
            degree.index_add_(0, dst, torch.ones_like(dst, dtype=x.dtype))
            degree = degree.clamp(min=1.0).unsqueeze(1)
            mean_neighbors = aggregated / degree
            return self.self_linear(x) + self.neighbor_linear(mean_neighbors)

    class SimpleGraphSAGE(_GraphRuntimeMixin, nn.Module):
        backend_name = "torch"
        model_name = "graphsage"

        def __init__(self) -> None:
            super().__init__()
            if num_layers < 1:
                raise ValueError("num_layers must be >= 1")
            self.layers = nn.ModuleList()
            dims = [input_dim] + [hidden_dim] * max(num_layers - 1, 0) + [output_dim]
            for in_dim, out_dim in zip(dims[:-1], dims[1:]):
                self.layers.append(GraphSAGELayer(in_dim, out_dim))
            self.activation = nn.ReLU()
            self.dropout = nn.Dropout(dropout)

        def forward(self, x, edge_index):
            hidden = x
            for layer_index, layer in enumerate(self.layers):
                hidden = layer(hidden, edge_index)
                if layer_index < len(self.layers) - 1:
                    hidden = self.activation(hidden)
                    hidden = self.dropout(hidden)
            return hidden

        def reconstruction_loss(self, embeddings, positive_edge_index):
            from torch_geometric.utils import negative_sampling

            pos_src, pos_dst = positive_edge_index
            pos_score = (embeddings[pos_src] * embeddings[pos_dst]).sum(dim=1)
            neg_edge_index = negative_sampling(
                positive_edge_index,
                num_nodes=embeddings.shape[0],
                num_neg_samples=positive_edge_index.shape[1],
                method="sparse",
            )
            neg_src, neg_dst = neg_edge_index
            neg_score = (embeddings[neg_src] * embeddings[neg_dst]).sum(dim=1)
            loss_pos = nn.functional.binary_cross_entropy_with_logits(pos_score, torch.ones_like(pos_score))
            loss_neg = nn.functional.binary_cross_entropy_with_logits(neg_score, torch.zeros_like(neg_score))
            return loss_pos + loss_neg

    return SimpleGraphSAGE()


def _build_pyg_graphsage_model(
    input_dim: int,
    hidden_dim: int = 64,
    output_dim: int = 64,
    num_layers: int = 2,
    dropout: float = 0.2,
):
    torch, nn, F = _require_torch()
    from torch_geometric.nn import SAGEConv
    from torch_geometric.utils import negative_sampling

    class PyGGraphSAGE(_GraphRuntimeMixin, nn.Module):
        backend_name = "pyg"
        model_name = "graphsage"

        def __init__(self) -> None:
            super().__init__()
            if num_layers < 1:
                raise ValueError("num_layers must be >= 1")
            self.convs = nn.ModuleList()
            dims = [input_dim] + [hidden_dim] * max(num_layers - 1, 0) + [output_dim]
            for in_dim, out_dim in zip(dims[:-1], dims[1:]):
                self.convs.append(SAGEConv(in_dim, out_dim))
            self.dropout = dropout

        def forward(self, x, edge_index):
            hidden = x
            for layer_index, conv in enumerate(self.convs):
                hidden = conv(hidden, edge_index)
                if layer_index < len(self.convs) - 1:
                    hidden = F.relu(hidden)
                    hidden = F.dropout(hidden, p=self.dropout, training=self.training)
            return hidden

        def reconstruction_loss(self, embeddings, positive_edge_index):
            pos_src, pos_dst = positive_edge_index
            pos_score = (embeddings[pos_src] * embeddings[pos_dst]).sum(dim=1)
            neg_edge_index = negative_sampling(
                positive_edge_index,
                num_nodes=embeddings.shape[0],
                num_neg_samples=positive_edge_index.shape[1],
                method="sparse",
            )
            neg_src, neg_dst = neg_edge_index
            neg_score = (embeddings[neg_src] * embeddings[neg_dst]).sum(dim=1)
            loss_pos = F.binary_cross_entropy_with_logits(pos_score, torch.ones_like(pos_score))
            loss_neg = F.binary_cross_entropy_with_logits(neg_score, torch.zeros_like(neg_score))
            return loss_pos + loss_neg

    return PyGGraphSAGE()


def build_graphsage_model(
    input_dim: int,
    hidden_dim: int = 64,
    output_dim: int = 64,
    num_layers: int = 2,
    dropout: float = 0.2,
    backend: str = "auto",
):
    status = inspect_graph_dependencies()
    normalized_backend = str(backend).strip().lower()

    if normalized_backend in {"auto", "pyg"} and status.torch_geometric_available:
        return _build_pyg_graphsage_model(
            input_dim=input_dim,
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            num_layers=num_layers,
            dropout=dropout,
        )

    if normalized_backend in {"auto", "torch"}:
        return _build_torch_graphsage_model(
            input_dim=input_dim,
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            num_layers=num_layers,
            dropout=dropout,
        )

    raise RuntimeError(
        f"Unable to build GraphSAGE with backend='{backend}'. Dependency inspection: {'; '.join(status.notes)}"
    )


def build_gcn_model(
    input_dim: int,
    hidden_dim: int = 64,
    output_dim: int = 64,
    num_layers: int = 2,
    dropout: float = 0.2,
    backend: str = "auto",
):
    status = inspect_graph_dependencies()
    normalized_backend = str(backend).strip().lower()
    if normalized_backend in {"auto", "pyg"} and status.torch_geometric_available:
        torch, nn, F = _require_torch()
        from torch_geometric.nn import GCNConv
        from torch_geometric.utils import negative_sampling

        class PyGGCN(_GraphRuntimeMixin, nn.Module):
            backend_name = "pyg"
            model_name = "gcn"

            def __init__(self) -> None:
                super().__init__()
                self.convs = nn.ModuleList()
                dims = [input_dim] + [hidden_dim] * max(num_layers - 1, 0) + [output_dim]
                for in_dim, out_dim in zip(dims[:-1], dims[1:]):
                    self.convs.append(GCNConv(in_dim, out_dim))
                self.dropout = dropout

            def forward(self, x, edge_index):
                hidden = x
                for layer_index, conv in enumerate(self.convs):
                    hidden = conv(hidden, edge_index)
                    if layer_index < len(self.convs) - 1:
                        hidden = F.relu(hidden)
                        hidden = F.dropout(hidden, p=self.dropout, training=self.training)
                return hidden

            def reconstruction_loss(self, embeddings, positive_edge_index):
                pos_src, pos_dst = positive_edge_index
                pos_score = (embeddings[pos_src] * embeddings[pos_dst]).sum(dim=1)
                neg_edge_index = negative_sampling(
                    positive_edge_index,
                    num_nodes=embeddings.shape[0],
                    num_neg_samples=positive_edge_index.shape[1],
                    method="sparse",
                )
                neg_src, neg_dst = neg_edge_index
                neg_score = (embeddings[neg_src] * embeddings[neg_dst]).sum(dim=1)
                loss_pos = F.binary_cross_entropy_with_logits(pos_score, torch.ones_like(pos_score))
                loss_neg = F.binary_cross_entropy_with_logits(neg_score, torch.zeros_like(neg_score))
                return loss_pos + loss_neg

        return PyGGCN()

    raise RuntimeError(
        f"GCN prototype currently expects torch_geometric. Dependency inspection: {'; '.join(status.notes)}"
    )


def build_gat_model(
    input_dim: int,
    hidden_dim: int = 64,
    output_dim: int = 64,
    num_layers: int = 2,
    dropout: float = 0.2,
    backend: str = "auto",
    heads: int = 2,
):
    status = inspect_graph_dependencies()
    normalized_backend = str(backend).strip().lower()
    if normalized_backend in {"auto", "pyg"} and status.torch_geometric_available:
        torch, nn, F = _require_torch()
        from torch_geometric.nn import GATConv
        from torch_geometric.utils import negative_sampling

        class PyGGAT(_GraphRuntimeMixin, nn.Module):
            backend_name = "pyg"
            model_name = "gat"

            def __init__(self) -> None:
                super().__init__()
                self.convs = nn.ModuleList()
                current_in = input_dim
                for layer_index in range(max(num_layers - 1, 1)):
                    self.convs.append(GATConv(current_in, hidden_dim, heads=heads, dropout=dropout))
                    current_in = hidden_dim * heads
                self.output_conv = GATConv(current_in, output_dim, heads=1, concat=False, dropout=dropout)
                self.dropout = dropout

            def forward(self, x, edge_index):
                hidden = x
                for conv in self.convs:
                    hidden = conv(hidden, edge_index)
                    hidden = F.elu(hidden)
                    hidden = F.dropout(hidden, p=self.dropout, training=self.training)
                hidden = self.output_conv(hidden, edge_index)
                return hidden

            def reconstruction_loss(self, embeddings, positive_edge_index):
                pos_src, pos_dst = positive_edge_index
                pos_score = (embeddings[pos_src] * embeddings[pos_dst]).sum(dim=1)
                neg_edge_index = negative_sampling(
                    positive_edge_index,
                    num_nodes=embeddings.shape[0],
                    num_neg_samples=positive_edge_index.shape[1],
                    method="sparse",
                )
                neg_src, neg_dst = neg_edge_index
                neg_score = (embeddings[neg_src] * embeddings[neg_dst]).sum(dim=1)
                loss_pos = F.binary_cross_entropy_with_logits(pos_score, torch.ones_like(pos_score))
                loss_neg = F.binary_cross_entropy_with_logits(neg_score, torch.zeros_like(neg_score))
                return loss_pos + loss_neg

        return PyGGAT()

    raise RuntimeError(
        f"GAT prototype currently expects torch_geometric. Dependency inspection: {'; '.join(status.notes)}"
    )


def build_graph_model(model_name: str, **kwargs):
    normalized = str(model_name).strip().lower()
    if normalized == "graphsage":
        return build_graphsage_model(**kwargs)
    if normalized == "gcn":
        return build_gcn_model(**kwargs)
    if normalized == "gat":
        return build_gat_model(**kwargs)
    raise ValueError(f"Unsupported graph model: {model_name}")
