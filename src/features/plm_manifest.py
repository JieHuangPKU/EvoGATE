"""
Minimal PLM manifest helpers for Phase 2 slim.
"""

import os

import pandas as pd

from src.features.plm_loaders import get_plm_dim


def build_plm_manifest(plm_name, source_path, key_type, n_keys, sample_keys):
    return pd.DataFrame(
        [
            {
                "plm_name": plm_name,
                "source_path": source_path,
                "file_format": os.path.splitext(source_path)[1].lstrip("."),
                "embedding_dim": get_plm_dim(plm_name),
                "key_type": key_type,
                "num_embedding_keys": n_keys,
                "sample_keys": ",".join(sample_keys[:5]),
            }
        ]
    )
