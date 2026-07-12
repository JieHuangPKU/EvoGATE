import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F


def _load_node2vec_class():
    try:
        from torch_geometric.nn import Node2Vec

        return Node2Vec
    except ImportError:
        from torch_geometric.nn.models import Node2Vec

        return Node2Vec


def set_seed(seed):
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def normalize_edge_index(edge_array):
    if edge_array.ndim != 2:
        raise ValueError(f"Expected 2D edge array, got shape {edge_array.shape}")
    if edge_array.shape[1] != 2 and edge_array.shape[0] == 2:
        edge_array = edge_array.T
    if edge_array.shape[1] != 2:
        raise ValueError(f"Expected edge array with two columns, got shape {edge_array.shape}")
    edge_index = torch.as_tensor(edge_array.T.astype(np.int64, copy=False), dtype=torch.long)
    return edge_index.contiguous()


def _normalize_backend_name(backend):
    raw = str(backend or "auto").strip().lower()
    aliases = {
        "auto": "auto",
        "pyg": "pyg",
        "torch_geometric": "pyg",
        "native_walk": "native_walk",
        "native": "native_walk",
        "walk": "native_walk",
        "svd": "svd",
        "spectral": "svd",
        "truncated_svd": "svd",
        "legacy_safe": "svd",
    }
    if raw not in aliases:
        raise ValueError("Unsupported node2vec backend '{}'".format(backend))
    return aliases[raw]


def _normalize_graph_contract(graph_contract):
    normalized = str(graph_contract or "undirected_symmetrized").strip().lower()
    allowed = {"directed_raw", "undirected_symmetrized"}
    if normalized not in allowed:
        raise ValueError("Unsupported graph contract '{}'".format(graph_contract))
    return normalized


def _resolve_backend(backend, require_true_node2vec):
    normalized_backend = _normalize_backend_name(backend)
    if normalized_backend == "auto":
        return "native_walk" if bool(require_true_node2vec) else "svd"
    if normalized_backend == "svd" and bool(require_true_node2vec):
        raise RuntimeError("require_true_node2vec=true forbids backend=svd")
    return normalized_backend


def _train_with_pyg_node2vec(edge_array, num_nodes, params, seed, device="cpu"):
    set_seed(seed)
    Node2Vec = _load_node2vec_class()
    edge_index = normalize_edge_index(edge_array)
    device_obj = torch.device(device)
    model = Node2Vec(
        edge_index=edge_index,
        embedding_dim=int(params["embedding_dim"]),
        walk_length=int(params["walk_length"]),
        context_size=int(params["context_size"]),
        walks_per_node=int(params["walks_per_node"]),
        p=float(params["p"]),
        q=float(params["q"]),
        num_negative_samples=int(params["num_negative_samples"]),
        sparse=True,
        num_nodes=int(num_nodes),
    ).to(device_obj)
    loader = model.loader(batch_size=int(params["batch_size"]), shuffle=True, num_workers=0)
    optimizer = torch.optim.SparseAdam(model.parameters(), lr=float(params["learning_rate"]))

    history_rows = []
    for epoch in range(int(params["epochs"])):
        model.train()
        epoch_losses = []
        for pos_rw, neg_rw in loader:
            optimizer.zero_grad()
            loss = model.loss(pos_rw.to(device_obj), neg_rw.to(device_obj))
            loss.backward()
            optimizer.step()
            epoch_losses.append(float(loss.item()))
        history_rows.append(
            {
                "epoch": epoch + 1,
                "mean_loss": float(np.mean(epoch_losses)) if epoch_losses else float("nan"),
                "batches": int(len(epoch_losses)),
            }
        )

    model.eval()
    with torch.no_grad():
        embeddings = model().detach().cpu().numpy().astype(np.float32, copy=False)

    if embeddings.shape != (int(num_nodes), int(params["embedding_dim"])):
        raise RuntimeError(
            "Unexpected node2vec embedding shape {}; expected {}".format(
                embeddings.shape, (int(num_nodes), int(params["embedding_dim"]))
            )
        )
    history_df = pd.DataFrame(history_rows)
    history_df["backend"] = "pyg_node2vec"
    history_df["fallback_used"] = False
    return embeddings, history_df


class _SkipGramNode2Vec(nn.Module):
    def __init__(self, num_nodes, embedding_dim):
        super().__init__()
        self.target_embeddings = nn.Embedding(int(num_nodes), int(embedding_dim), sparse=True)
        self.context_embeddings = nn.Embedding(int(num_nodes), int(embedding_dim), sparse=True)
        nn.init.xavier_uniform_(self.target_embeddings.weight)
        nn.init.xavier_uniform_(self.context_embeddings.weight)

    def forward(self):
        return self.target_embeddings.weight


def _build_neighbors(edge_array, num_nodes):
    neighbors = [[] for _ in range(int(num_nodes))]
    edge_array = np.asarray(edge_array, dtype=np.int64)
    for src, dst in edge_array.tolist():
        if src < 0 or dst < 0 or src >= int(num_nodes) or dst >= int(num_nodes):
            raise ValueError("edge_index contains out-of-range node ids")
        neighbors[int(src)].append(int(dst))
    return [np.asarray(row, dtype=np.int64) for row in neighbors]


def _build_transition_weights(neighbors, p, q):
    adjacency_sets = [set(row.tolist()) for row in neighbors]
    edge_weights = {}
    inv_p = 1.0 / float(p)
    inv_q = 1.0 / float(q)
    for prev_node, row in enumerate(neighbors):
        for current_node in row.tolist():
            current_neighbors = neighbors[current_node]
            if current_neighbors.size == 0:
                continue
            weights = np.ones(current_neighbors.shape[0], dtype=np.float64)
            previous_adj = adjacency_sets[prev_node]
            for idx, candidate in enumerate(current_neighbors.tolist()):
                if candidate == prev_node:
                    weights[idx] = inv_p
                elif candidate not in previous_adj:
                    weights[idx] = inv_q
            weight_sum = float(weights.sum())
            if weight_sum <= 0.0:
                weights = np.full(current_neighbors.shape[0], 1.0 / float(current_neighbors.shape[0]), dtype=np.float64)
            else:
                weights /= weight_sum
            edge_weights[(prev_node, current_node)] = weights.astype(np.float64, copy=False)
    return edge_weights


def _sample_walk(start_nodes, neighbors, transition_weights, walk_length, rng):
    walks = np.full((int(start_nodes.shape[0]), int(walk_length) + 1), -1, dtype=np.int64)
    walks[:, 0] = start_nodes
    for row_idx, start_node in enumerate(start_nodes.tolist()):
        current = int(start_node)
        previous = -1
        for step_idx in range(1, int(walk_length) + 1):
            current_neighbors = neighbors[current]
            if current_neighbors.size == 0:
                walks[row_idx, step_idx:] = current
                break
            if previous < 0:
                choice_idx = int(rng.randint(0, current_neighbors.shape[0]))
            else:
                probabilities = transition_weights.get((previous, current))
                if probabilities is None:
                    choice_idx = int(rng.randint(0, current_neighbors.shape[0]))
                else:
                    choice_idx = int(rng.choice(current_neighbors.shape[0], p=probabilities))
            nxt = int(current_neighbors[choice_idx])
            walks[row_idx, step_idx] = nxt
            previous, current = current, nxt
    return walks


def _iter_positive_pairs(walks, context_size, rng):
    walk_width = int(walks.shape[1])
    window = max(1, int(context_size))
    for walk in walks:
        for center_idx in range(walk_width):
            center = int(walk[center_idx])
            if center < 0:
                continue
            left = max(0, center_idx - window)
            right = min(walk_width, center_idx + window + 1)
            candidates = [idx for idx in range(left, right) if idx != center_idx and int(walk[idx]) >= 0]
            if not candidates:
                continue
            context_idx = int(candidates[int(rng.randint(0, len(candidates)))])
            yield center, int(walk[context_idx])


def _train_with_native_walk_node2vec(edge_array, num_nodes, params, seed, device="cpu"):
    set_seed(seed)
    device_obj = torch.device(device)
    num_nodes = int(num_nodes)
    if num_nodes <= 0:
        raise ValueError("num_nodes must be positive")
    if float(params["p"]) <= 0.0 or float(params["q"]) <= 0.0:
        raise ValueError("Node2Vec p and q must both be > 0")
    edge_array = np.asarray(edge_array, dtype=np.int64)
    if edge_array.size == 0:
        raise RuntimeError("true walk-based Node2Vec requires a non-empty graph")

    neighbors = _build_neighbors(edge_array, num_nodes)
    if all(row.size == 0 for row in neighbors):
        raise RuntimeError("true walk-based Node2Vec requires at least one reachable edge")
    transition_weights = _build_transition_weights(neighbors, float(params["p"]), float(params["q"]))
    model = _SkipGramNode2Vec(num_nodes=num_nodes, embedding_dim=int(params["embedding_dim"])).to(device_obj)
    optimizer = torch.optim.SparseAdam(model.parameters(), lr=float(params["learning_rate"]))
    node_ids = np.arange(num_nodes, dtype=np.int64)
    negative_count = max(1, int(params["num_negative_samples"]))
    batch_size = max(1, int(params["batch_size"]))
    pair_batch_size = max(batch_size, batch_size * max(1, int(params["walk_length"])))
    rng = np.random.RandomState(int(seed))
    history_rows = []

    for epoch in range(1, int(params["epochs"]) + 1):
        epoch_rng = np.random.RandomState(int(seed) + epoch)
        start_nodes = np.repeat(node_ids, int(params["walks_per_node"]))
        epoch_rng.shuffle(start_nodes)
        epoch_losses = []
        batch_pairs = []
        for offset in range(0, start_nodes.shape[0], batch_size):
            walk_batch = _sample_walk(
                start_nodes=start_nodes[offset:offset + batch_size],
                neighbors=neighbors,
                transition_weights=transition_weights,
                walk_length=int(params["walk_length"]),
                rng=epoch_rng,
            )
            for pair in _iter_positive_pairs(walk_batch, int(params["context_size"]), epoch_rng):
                batch_pairs.append(pair)
                if len(batch_pairs) >= pair_batch_size:
                    centers = torch.as_tensor([pair_item[0] for pair_item in batch_pairs], dtype=torch.long, device=device_obj)
                    contexts = torch.as_tensor([pair_item[1] for pair_item in batch_pairs], dtype=torch.long, device=device_obj)
                    negatives = torch.as_tensor(
                        rng.randint(0, num_nodes, size=(len(batch_pairs), negative_count)),
                        dtype=torch.long,
                        device=device_obj,
                    )
                    optimizer.zero_grad()
                    center_vec = model.target_embeddings(centers)
                    pos_vec = model.context_embeddings(contexts)
                    neg_vec = model.context_embeddings(negatives)
                    pos_score = (center_vec * pos_vec).sum(dim=1)
                    neg_score = torch.bmm(neg_vec, center_vec.unsqueeze(-1)).squeeze(-1)
                    loss = -F.logsigmoid(pos_score).mean() - F.logsigmoid(-neg_score).mean()
                    loss.backward()
                    optimizer.step()
                    epoch_losses.append(float(loss.item()))
                    batch_pairs = []
        if batch_pairs:
            centers = torch.as_tensor([pair_item[0] for pair_item in batch_pairs], dtype=torch.long, device=device_obj)
            contexts = torch.as_tensor([pair_item[1] for pair_item in batch_pairs], dtype=torch.long, device=device_obj)
            negatives = torch.as_tensor(
                rng.randint(0, num_nodes, size=(len(batch_pairs), negative_count)),
                dtype=torch.long,
                device=device_obj,
            )
            optimizer.zero_grad()
            center_vec = model.target_embeddings(centers)
            pos_vec = model.context_embeddings(contexts)
            neg_vec = model.context_embeddings(negatives)
            pos_score = (center_vec * pos_vec).sum(dim=1)
            neg_score = torch.bmm(neg_vec, center_vec.unsqueeze(-1)).squeeze(-1)
            loss = -F.logsigmoid(pos_score).mean() - F.logsigmoid(-neg_score).mean()
            loss.backward()
            optimizer.step()
            epoch_losses.append(float(loss.item()))
            batch_pairs = []
        history_rows.append(
            {
                "epoch": epoch,
                "mean_loss": float(np.mean(epoch_losses)) if epoch_losses else float("nan"),
                "batches": int(len(epoch_losses)),
                "backend": "native_walk_node2vec",
                "fallback_used": False,
            }
        )

    with torch.no_grad():
        embeddings = model().detach().cpu().numpy().astype(np.float32, copy=False)
    if embeddings.shape != (num_nodes, int(params["embedding_dim"])):
        raise RuntimeError("Unexpected native walk node2vec embedding shape {}".format(embeddings.shape))
    return embeddings, pd.DataFrame(history_rows)


def _train_with_svd_fallback(edge_array, num_nodes, params, seed):
    raise RuntimeError(
        "SVD fallback is disabled for the clean true-Node2Vec benchmark path. "
        "Use backend=native_walk or backend=pyg, or set require_true_node2vec=false outside benchmark mode."
    )


def train_node2vec_embeddings(
    edge_array,
    num_nodes,
    params,
    seed,
    device="cpu",
    backend="auto",
    require_true_node2vec=False,
):
    resolved_backend = _resolve_backend(backend, require_true_node2vec=require_true_node2vec)
    if resolved_backend == "svd":
        return _train_with_svd_fallback(edge_array, num_nodes, params, seed)
    if resolved_backend == "pyg":
        try:
            return _train_with_pyg_node2vec(edge_array, num_nodes, params, seed, device=device)
        except Exception as exc:
            raise RuntimeError("PyG Node2Vec failed under require_true_node2vec={}: {}".format(require_true_node2vec, exc))
    if resolved_backend == "native_walk":
        return _train_with_native_walk_node2vec(edge_array, num_nodes, params, seed, device=device)
    raise ValueError("Unsupported resolved backend '{}'".format(resolved_backend))


def write_node2vec_outputs(embeddings, params, output_dir, history_df, metadata=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = dict(metadata or {})
    backend_name = str(metadata.get("embedding_backend", "unknown"))
    fallback_used = bool(metadata.get("fallback_used", False))
    graph_contract = _normalize_graph_contract(metadata.get("graph_contract", "undirected_symmetrized"))
    np.save(output_dir / "node2vec_embeddings.npy", embeddings)
    pd.DataFrame(
        [
            {
                "feature_block": "node2vec",
                "start_col": 0,
                "end_col": int(embeddings.shape[1] - 1),
                "dimension": int(embeddings.shape[1]),
                "data_source": "ppi_node2vec",
                "missing_strategy": "not_applicable",
            }
        ]
    ).to_csv(output_dir / "feature_schema.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "embedding_method": str(metadata.get("embedding_method", "node2vec_walk_based")),
                "embedding_backend": backend_name,
                "fallback_used": str(fallback_used).lower(),
                "feature_block": "node2vec",
                "feature_dim": int(embeddings.shape[1]),
                "dimensions": int(params["embedding_dim"]),
                "embedding_dim": int(params["embedding_dim"]),
                "walk_length": int(params["walk_length"]),
                "context_size": int(params["context_size"]),
                "num_walks": int(params["walks_per_node"]),
                "walks_per_node": int(params["walks_per_node"]),
                "num_negative_samples": int(params["num_negative_samples"]),
                "epochs": int(params["epochs"]),
                "batch_size": int(params["batch_size"]),
                "learning_rate": float(params["learning_rate"]),
                "p": float(params["p"]),
                "q": float(params["q"]),
                "graph_contract": graph_contract,
                "require_true_node2vec": str(bool(metadata.get("require_true_node2vec", False))).lower(),
            }
        ]
    ).to_csv(output_dir / "node2vec_summary.tsv", sep="\t", index=False)
    history_df.to_csv(output_dir / "node2vec_training_log.tsv", sep="\t", index=False)


def parse_args():
    parser = argparse.ArgumentParser(description="Train node2vec embeddings on a benchmark graph")
    parser.add_argument("--edge-index-path", required=True, type=str)
    parser.add_argument("--num-nodes", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=str)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--embedding-dim", default=64, type=int)
    parser.add_argument("--walk-length", default=32, type=int)
    parser.add_argument("--context-size", default=32, type=int)
    parser.add_argument("--walks-per-node", default=16, type=int)
    parser.add_argument("--num-negative-samples", default=1, type=int)
    parser.add_argument("--epochs", default=50, type=int)
    parser.add_argument("--batch-size", default=256, type=int)
    parser.add_argument("--learning-rate", default=0.005, type=float)
    parser.add_argument("--p", default=1.0, type=float)
    parser.add_argument("--q", default=1.0, type=float)
    parser.add_argument("--device", default="cpu", type=str)
    parser.add_argument("--backend", default="auto", type=str)
    parser.add_argument("--require-true-node2vec", action="store_true")
    parser.add_argument("--graph-contract", default="undirected_symmetrized", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    params = {
        "embedding_dim": args.embedding_dim,
        "walk_length": args.walk_length,
        "context_size": args.context_size,
        "walks_per_node": args.walks_per_node,
        "num_negative_samples": args.num_negative_samples,
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "learning_rate": args.learning_rate,
        "p": args.p,
        "q": args.q,
    }
    edge_array = np.load(args.edge_index_path, allow_pickle=True)
    embeddings, history_df = train_node2vec_embeddings(
        edge_array=edge_array,
        num_nodes=args.num_nodes,
        params=params,
        seed=int(args.seed),
        device=args.device,
        backend=args.backend,
        require_true_node2vec=bool(args.require_true_node2vec),
    )
    backend_name = "unknown"
    if not history_df.empty and "backend" in history_df.columns:
        backend_name = str(history_df["backend"].iloc[0])
    write_node2vec_outputs(
        embeddings,
        params,
        args.output_dir,
        history_df,
        metadata={
            "embedding_method": "node2vec_walk_based",
            "embedding_backend": backend_name,
            "fallback_used": False,
            "require_true_node2vec": bool(args.require_true_node2vec),
            "graph_contract": args.graph_contract,
        },
    )


if __name__ == "__main__":
    main()
