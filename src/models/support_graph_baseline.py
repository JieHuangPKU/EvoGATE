from src.models.gnn_graphsage import SupportGraphSAGE
from src.models.gnn_gcn import SupportGCN
from src.models.gnn_gat import SupportGAT


def build_support_graph_model(model_type, in_feats, hidden_feats, out_feats, num_layers=2, dropout=0.2, num_heads=2):
    model_type = str(model_type).strip().lower()
    if model_type == "graphsage":
        return SupportGraphSAGE(in_feats, hidden_feats, out_feats, num_layers, dropout)
    if model_type == "gcn":
        return SupportGCN(in_feats, hidden_feats, out_feats, num_layers, dropout)
    if model_type == "gat":
        return SupportGAT(in_feats, hidden_feats, out_feats, num_layers, dropout, num_heads)
    raise ValueError("错误：不支持的图模型类型 {}".format(model_type))
