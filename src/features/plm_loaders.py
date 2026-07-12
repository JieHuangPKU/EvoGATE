"""
Minimal Phase 2 slim PLM loaders for ESM2 and ProtT5.
"""

import os

import h5py
import numpy as np


SUPPORTED_PLM = {
    "esm2": {"dim": 1280},
    "prott5": {"dim": 1024},
}


def load_h5_embeddings(path):
    if not os.path.exists(path):
        raise FileNotFoundError("Missing embedding file: {}".format(path))
    with h5py.File(path, "r") as handle:
        return {key: handle[key][:] for key in handle.keys()}


def get_plm_dim(plm_name):
    if plm_name not in SUPPORTED_PLM:
        raise ValueError("Unsupported PLM: {}".format(plm_name))
    return int(SUPPORTED_PLM[plm_name]["dim"])


def zero_vector(plm_name):
    return np.zeros((get_plm_dim(plm_name),), dtype=np.float32)
