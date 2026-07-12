"""
Original EPGAT-compatible GAT model.

Migration source:
- /home/jiehuang/software/fungi/EPGAT/models/gat/gat_pytorch.py

Boundary:
- This file is an original-compatible adapter for Phase 1.
- It is not the unified upgraded ProGATE_v2 graph model layer.
"""

import warnings

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Parameter
from torch_geometric.nn.conv import MessagePassing
from torch_geometric.nn.inits import glorot, zeros
from torch_geometric.utils import add_self_loops, remove_self_loops, softmax


class EPGATOriginalConv(MessagePassing):
    def __init__(self, in_channels, out_channels, heads=1, concat=True,
                 negative_slope=0.2, dropout=0.0, bias=True, **kwargs):
        super(EPGATOriginalConv, self).__init__(aggr="add", **kwargs)
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.heads = heads
        self.concat = concat
        self.negative_slope = negative_slope
        self.dropout = dropout
        self.weight = Parameter(torch.Tensor(in_channels, heads * out_channels))
        self.att = Parameter(torch.Tensor(1, heads, 2 * out_channels))
        if bias and concat:
            self.bias = Parameter(torch.Tensor(heads * out_channels))
        elif bias and not concat:
            self.bias = Parameter(torch.Tensor(out_channels))
        else:
            self.register_parameter("bias", None)
        self.reset_parameters()

    def reset_parameters(self):
        glorot(self.weight)
        glorot(self.att)
        zeros(self.bias)

    def forward(self, x, edge_index, size=None, edge_weights=None, return_alpha=False):
        self.return_alpha = return_alpha
        self.edge_weights = edge_weights
        if size is None and torch.is_tensor(x):
            edge_index, _ = remove_self_loops(edge_index)
            edge_index, _ = add_self_loops(edge_index, num_nodes=x.size(self.node_dim))
        if torch.is_tensor(x):
            x = torch.matmul(x, self.weight)
        else:
            x = (
                None if x[0] is None else torch.matmul(x[0], self.weight),
                None if x[1] is None else torch.matmul(x[1], self.weight),
            )
        if self.edge_weights is not None:
            self.edge_weights = self.edge_weights.view(-1)
        if self.return_alpha:
            return self.propagate(edge_index, size=size, x=x, edge_weight=self.edge_weights), self.alpha, edge_index
        return self.propagate(edge_index, size=size, x=x, edge_weight=self.edge_weights)

    def message(self, edge_index_i, x_i, x_j, size_i):
        x_j = x_j.view(-1, self.heads, self.out_channels)
        if x_i is None:
            alpha = (x_j * self.att[:, :, self.out_channels:]).sum(dim=-1)
        else:
            x_i = x_i.view(-1, self.heads, self.out_channels)
            alpha = (torch.cat([x_i, x_j], dim=-1) * self.att).sum(dim=-1)
        alpha = F.leaky_relu(alpha, self.negative_slope)
        alpha = softmax(alpha, edge_index_i, num_nodes=size_i)
        if self.return_alpha:
            self.alpha = alpha
        alpha = F.dropout(alpha, p=self.dropout, training=self.training)
        if self.edge_weights is not None:
            edge_weights = self.edge_weights.reshape((-1, 1, 1))
            return x_j * alpha.view(-1, self.heads, 1) * edge_weights
        return x_j * alpha.view(-1, self.heads, 1)

    def update(self, aggr_out):
        if self.concat:
            aggr_out = aggr_out.view(-1, self.heads * self.out_channels)
        else:
            aggr_out = aggr_out.mean(dim=1)
        if self.bias is not None:
            aggr_out = aggr_out + self.bias
        return aggr_out


class EPGATOriginalGAT(nn.Module):
    def __init__(self, in_feats=1, h_feats=None, heads=None, dropout=0.6, negative_slope=0.2, linear_layer=None):
        super(EPGATOriginalGAT, self).__init__()
        h_feats = h_feats or [8, 1]
        heads = heads or [8, 1]
        self.dropout = dropout
        self.linear_layer = linear_layer
        self.layers = nn.ModuleList()
        if self.linear_layer is not None:
            self.linear = nn.Linear(in_feats, linear_layer)
        in_size = in_feats if linear_layer is None else linear_layer
        for idx, h_feat in enumerate(h_feats):
            last = idx + 1 == len(h_feats)
            self.layers.append(
                EPGATOriginalConv(
                    in_size,
                    h_feat,
                    heads=heads[idx],
                    dropout=dropout,
                    negative_slope=negative_slope,
                    concat=False if last else True,
                )
            )
            in_size = h_feat * heads[idx]

    def forward(self, x, edge_index, edge_weights=None, return_alphas=False):
        if self.linear_layer is not None:
            x = self.linear(x)
        alphas = []
        for idx, layer in enumerate(self.layers):
            try:
                if idx < len(self.layers) - 1:
                    if return_alphas:
                        x, alpha, _ = layer(x, edge_index, edge_weights=edge_weights, return_alpha=True)
                        alphas.append(alpha)
                    else:
                        x = layer(x, edge_index, edge_weights=edge_weights)
                    x = F.relu(x)
                    x = F.dropout(x, self.dropout, training=self.training)
                else:
                    if return_alphas:
                        x, alpha, edge_idx = layer(x, edge_index, edge_weights=edge_weights, return_alpha=True)
                        alphas.append(alpha)
                        if x.dim() == 1:
                            x = x.unsqueeze(-1)
                        return x, alphas, edge_idx
                    else:
                        x = layer(x, edge_index, edge_weights=edge_weights)
                        if x.dim() == 1:
                            x = x.unsqueeze(-1)
            except RuntimeError as exc:
                warnings.warn("Layer {} failed: {}".format(idx, exc))
                raise
        return x
