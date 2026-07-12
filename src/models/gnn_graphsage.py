import torch.nn as nn
import torch.nn.functional as F

try:
    from dgl.nn.pytorch import SAGEConv
except ImportError:
    SAGEConv = None


class SupportGraphSAGE(nn.Module):
    def __init__(self, in_feats, hidden_feats, out_feats, num_layers, dropout):
        super(SupportGraphSAGE, self).__init__()
        if SAGEConv is None:
            raise ImportError("错误：当前环境缺少 dgl，请先执行 conda activate EPGAT。")
        self.layers = nn.ModuleList()
        self.layers.append(SAGEConv(in_feats, hidden_feats, "mean"))
        for _ in range(max(num_layers - 2, 0)):
            self.layers.append(SAGEConv(hidden_feats, hidden_feats, "mean"))
        self.layers.append(SAGEConv(hidden_feats, out_feats, "mean"))
        self.dropout = nn.Dropout(dropout)

    def forward(self, graph, features):
        h = features
        for idx, layer in enumerate(self.layers):
            h = layer(graph, h)
            if idx < len(self.layers) - 1:
                h = F.relu(h)
                h = self.dropout(h)
        return h
