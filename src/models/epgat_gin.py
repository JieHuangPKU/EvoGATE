"""
Original-compatible GIN for the archived replay benchmark.

Migration source:
- /home/jiehuang/software/fungi/EPGAT/models/gin/gin_pytorch.py
"""

import torch.nn as nn
from torch.nn import BatchNorm1d, Dropout, Linear, ReLU, Sequential
from torch_geometric.nn import GINConv
from torch_geometric.nn.inits import reset


class EPGATOriginalGIN(nn.Module):
    def __init__(self, in_feats, dim_h=64, dropout=0.5):
        super().__init__()
        self.dropout_layer = Dropout(dropout)

        nn1 = Sequential(Linear(in_feats, dim_h), BatchNorm1d(dim_h), ReLU(), Linear(dim_h, dim_h), BatchNorm1d(dim_h), ReLU())
        nn2 = Sequential(Linear(dim_h, dim_h), BatchNorm1d(dim_h), ReLU(), Linear(dim_h, dim_h), BatchNorm1d(dim_h), ReLU())
        nn3 = Sequential(Linear(dim_h, dim_h), BatchNorm1d(dim_h), ReLU(), Linear(dim_h, dim_h), BatchNorm1d(dim_h), ReLU())
        self.conv1 = GINConv(nn1)
        self.conv2 = GINConv(nn2)
        self.conv3 = GINConv(nn3)
        self.classifier = Linear(dim_h, 1)
        self.reset_parameters()

    def reset_parameters(self):
        reset(self.conv1)
        reset(self.conv2)
        reset(self.conv3)
        self.classifier.reset_parameters()

    def forward(self, x, edge_index, edge_weights=None):
        x = self.conv1(x, edge_index)
        x = self.dropout_layer(x)
        x = self.conv2(x, edge_index)
        x = self.dropout_layer(x)
        x = self.conv3(x, edge_index)
        x = self.dropout_layer(x)
        x = self.classifier(x)
        return x
