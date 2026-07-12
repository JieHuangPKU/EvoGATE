from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC, SVC


SUPPORTED_BASELINE_MODELS = {
    "logistic_regression",
    "mlp",
    "random_forest",
    "svm",
}

MODEL_ALIASES = {
    "logreg": "logistic_regression",
    "logistic": "logistic_regression",
    "logistic_regression": "logistic_regression",
    "mlp": "mlp",
    "rf": "random_forest",
    "random_forest": "random_forest",
    "svm": "svm",
}


def normalize_model_name(model_name: str) -> str:
    normalized = MODEL_ALIASES.get(str(model_name).strip().lower())
    if normalized is None:
        raise ValueError(
            f"Unsupported baseline model '{model_name}'. "
            f"Supported: {sorted(SUPPORTED_BASELINE_MODELS)}"
        )
    return normalized


@dataclass
class BaselineModelWrapper:
    model_name: str
    model: object
    model_params: dict[str, Any]

    def fit(self, x: np.ndarray, y: np.ndarray) -> None:
        self.model.fit(x, y)

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        if hasattr(self.model, "predict_proba"):
            probabilities = self.model.predict_proba(x)
            if probabilities.ndim != 2 or probabilities.shape[1] < 2:
                raise ValueError("Expected binary predict_proba output with two probability columns")
            return probabilities[:, 1]
        if hasattr(self.model, "decision_function"):
            scores = np.asarray(self.model.decision_function(x), dtype=np.float64)
            return 1.0 / (1.0 + np.exp(-scores))
        raise ValueError(f"Model '{self.model_name}' does not support probability inference")

    def predict(self, x: np.ndarray, threshold: float = 0.5) -> np.ndarray:
        return (self.predict_proba(x) >= threshold).astype(int)

    def save(self, output_path: str | Path) -> None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("wb") as handle:
            pickle.dump(
                {
                    "model_name": self.model_name,
                    "model": self.model,
                    "model_params": self.model_params,
                },
                handle,
            )

    @classmethod
    def load(cls, input_path: str | Path) -> "BaselineModelWrapper":
        with Path(input_path).open("rb") as handle:
            payload = pickle.load(handle)
        model_name = payload.get("model_name", payload.get("model_type", "logistic_regression"))
        return cls(
            model_name=normalize_model_name(model_name),
            model=payload["model"],
            model_params=dict(payload.get("model_params", {})),
        )


def _pipeline_with_scaler(model: object) -> Pipeline:
    return Pipeline(
        steps=[
            ("scaler", StandardScaler()),
            ("classifier", model),
        ]
    )


def create_baseline_model(
    model_name: str,
    config: dict[str, Any] | None = None,
    random_state: int = 42,
) -> BaselineModelWrapper:
    config = dict(config or {})
    model_name = normalize_model_name(model_name)

    if model_name == "logistic_regression":
        params = {
            "max_iter": int(config.get("max_iter", 2000)),
            "class_weight": config.get("class_weight", "balanced"),
            "solver": config.get("solver", "lbfgs"),
            "c": float(config.get("c", 1.0)),
        }
        model = _pipeline_with_scaler(
            LogisticRegression(
                max_iter=params["max_iter"],
                class_weight=params["class_weight"],
                solver=params["solver"],
                C=params["c"],
                random_state=random_state,
            )
        )
        return BaselineModelWrapper(model_name=model_name, model=model, model_params=params)

    if model_name == "mlp":
        raw_hidden = config.get("hidden_layer_sizes", [256, 64])
        if isinstance(raw_hidden, (list, tuple)):
            hidden_layer_sizes = tuple(int(value) for value in raw_hidden)
        else:
            hidden_layer_sizes = tuple(int(value) for value in str(raw_hidden).split(",") if str(value).strip())
        params = {
            "hidden_layer_sizes": hidden_layer_sizes,
            "activation": config.get("activation", "relu"),
            "solver": config.get("solver", "adam"),
            "batch_size": int(config.get("batch_size", 64)),
            "max_iter": int(config.get("max_iter", 300)),
            "learning_rate_init": float(config.get("learning_rate_init", 1e-3)),
            "alpha": float(config.get("alpha", 1e-4)),
            "early_stopping": bool(config.get("early_stopping", True)),
        }
        model = _pipeline_with_scaler(
            MLPClassifier(
                hidden_layer_sizes=params["hidden_layer_sizes"],
                activation=params["activation"],
                solver=params["solver"],
                batch_size=params["batch_size"],
                max_iter=params["max_iter"],
                learning_rate_init=params["learning_rate_init"],
                alpha=params["alpha"],
                early_stopping=params["early_stopping"],
                random_state=random_state,
            )
        )
        return BaselineModelWrapper(model_name=model_name, model=model, model_params=params)

    if model_name == "random_forest":
        params = {
            "n_estimators": int(config.get("n_estimators", 400)),
            "max_depth": None if config.get("max_depth", None) in [None, "", "null"] else int(config["max_depth"]),
            "min_samples_split": int(config.get("min_samples_split", 2)),
            "min_samples_leaf": int(config.get("min_samples_leaf", 1)),
            "class_weight": config.get("class_weight", "balanced"),
            "n_jobs": int(config.get("n_jobs", -1)),
        }
        model = RandomForestClassifier(
            n_estimators=params["n_estimators"],
            max_depth=params["max_depth"],
            min_samples_split=params["min_samples_split"],
            min_samples_leaf=params["min_samples_leaf"],
            class_weight=params["class_weight"],
            n_jobs=params["n_jobs"],
            random_state=random_state,
        )
        return BaselineModelWrapper(model_name=model_name, model=model, model_params=params)

    if model_name == "svm":
        params = {
            "c": float(config.get("c", 1.0)),
            "kernel": config.get("kernel", "rbf"),
            "gamma": config.get("gamma", "scale"),
            "class_weight": config.get("class_weight", "balanced"),
            "max_iter": int(config.get("max_iter", 5000)),
        }
        if params["kernel"] == "linear":
            classifier = LinearSVC(
                C=params["c"],
                class_weight=params["class_weight"],
                max_iter=params["max_iter"],
                random_state=random_state,
            )
        else:
            classifier = SVC(
                C=params["c"],
                kernel=params["kernel"],
                gamma=params["gamma"],
                class_weight=params["class_weight"],
                probability=False,
                random_state=random_state,
            )
        model = _pipeline_with_scaler(classifier)
        return BaselineModelWrapper(model_name=model_name, model=model, model_params=params)

    raise ValueError(f"Unsupported baseline model: {model_name}")
