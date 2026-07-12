"""
Original-compatible GCN for the archived replay benchmark.

Migration source:
- /home/jiehuang/software/fungi/EPGAT/models/gcn/gcn_pytorch.py
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv


class EPGATOriginalGCN(nn.Module):
    def __init__(self, in_feats, h_layers=None, dropout=0.5):
        super().__init__()
        h_layers = h_layers or [16, 1]
        self.dropout = dropout
        self.layers = nn.ModuleList()
        hidden_in = in_feats
        for out_feats in h_layers:
            self.layers.append(GCNConv(hidden_in, out_feats))
            hidden_in = out_feats

    def forward(self, x, edge_index, edge_weights=None):
        for layer in self.layers[:-1]:
            x = F.relu(layer(x, edge_index, edge_weight=edge_weights))
            x = F.dropout(x, self.dropout, training=self.training)
        x = self.layers[-1](x, edge_index, edge_weight=edge_weights)
        return x
