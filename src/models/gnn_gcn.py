import torch.nn as nn
import torch.nn.functional as F

try:
    from dgl.nn.pytorch import GraphConv
except ImportError:
    GraphConv = None


class SupportGCN(nn.Module):
    def __init__(self, in_feats, hidden_feats, out_feats, num_layers, dropout):
        super(SupportGCN, self).__init__()
        if GraphConv is None:
            raise ImportError("错误：当前环境缺少 dgl，请先执行 conda activate EPGAT。")
        self.layers = nn.ModuleList()
        self.layers.append(GraphConv(in_feats, hidden_feats))
        for _ in range(max(num_layers - 2, 0)):
            self.layers.append(GraphConv(hidden_feats, hidden_feats))
        self.layers.append(GraphConv(hidden_feats, out_feats))
        self.dropout = nn.Dropout(dropout)

    def forward(self, graph, features, edge_weight=None):
        h = features
        for idx, layer in enumerate(self.layers):
            try:
                h = layer(graph, h, edge_weight=edge_weight)
            except TypeError:
                h = layer(graph, h)
            if idx < len(self.layers) - 1:
                h = F.relu(h)
                h = self.dropout(h)
        return h
