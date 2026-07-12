import torch.nn as nn
import torch.nn.functional as F

try:
    from dgl.nn.pytorch import GATConv
except ImportError:
    GATConv = None


class SupportGAT(nn.Module):
    def __init__(self, in_feats, hidden_feats, out_feats, num_layers, dropout, num_heads):
        super(SupportGAT, self).__init__()
        if GATConv is None:
            raise ImportError("错误：当前环境缺少 dgl，请先执行 conda activate EPGAT。")
        self.layers = nn.ModuleList()
        self.layers.append(GATConv(in_feats, hidden_feats, num_heads, feat_drop=dropout, attn_drop=dropout))
        current_dim = hidden_feats * num_heads
        for _ in range(max(num_layers - 2, 0)):
            self.layers.append(GATConv(current_dim, hidden_feats, num_heads, feat_drop=dropout, attn_drop=dropout))
            current_dim = hidden_feats * num_heads
        self.out_layer = GATConv(current_dim, out_feats, 1, feat_drop=dropout, attn_drop=dropout)
        self.dropout = nn.Dropout(dropout)

    def forward(self, graph, features):
        h = features
        for layer in self.layers:
            h = layer(graph, h).flatten(1)
            h = F.elu(h)
            h = self.dropout(h)
        return self.out_layer(graph, h).mean(1)
