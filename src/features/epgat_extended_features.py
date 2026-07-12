"""
Phase 2 slim feature assembly helpers.
"""

import numpy as np
import pandas as pd

from src.features.epgat_legacy_features import zscore_matrix


def extend_feature_schema(base_schema, feature_block, start_col, dimension, data_source):
    end_col = start_col + dimension - 1
    extra = pd.DataFrame(
        [
            {
                "feature_block": feature_block,
                "start_col": start_col,
                "end_col": end_col,
                "dimension": dimension,
                "data_source": data_source,
                "missing_strategy": "key_map_fill_zero",
            }
        ]
    )
    return pd.concat([base_schema, extra], ignore_index=True)


def concatenate_feature_blocks(base_x, extra_x):
    merged = np.concatenate([base_x.astype(np.float32), extra_x.astype(np.float32)], axis=1)
    return zscore_matrix(merged).astype(np.float32)
