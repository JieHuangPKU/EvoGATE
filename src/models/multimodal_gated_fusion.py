import torch
import torch.nn as nn


class MultimodalGatedFusion(nn.Module):
    def __init__(
        self,
        omics_dim: int,
        esm2_dim: int,
        hidden_dim: int = 256,
        dropout: float = 0.2,
        fusion_mode: str = "gated",
    ):
        super().__init__()
        self.fusion_mode = str(fusion_mode).strip().lower()
        self.omics_dim = int(omics_dim)
        self.esm2_dim = int(esm2_dim)
        self.hidden_dim = int(hidden_dim)
        self.dropout = float(dropout)
        if self.fusion_mode not in {"gated", "residual_gated_concat"}:
            raise ValueError(f"Unsupported fusion_mode '{fusion_mode}'")

        if self.fusion_mode == "gated":
            self.omics_projection = nn.Sequential(
                nn.Linear(self.omics_dim, self.hidden_dim),
                nn.LayerNorm(self.hidden_dim),
                nn.GELU(),
                nn.Dropout(self.dropout),
            )
            self.esm2_projection = nn.Sequential(
                nn.Linear(self.esm2_dim, self.hidden_dim),
                nn.LayerNorm(self.hidden_dim),
                nn.GELU(),
                nn.Dropout(self.dropout),
            )
            self.gate_linear = nn.Linear(self.hidden_dim * 2, self.hidden_dim)
            self.output_dim = self.hidden_dim
        else:
            self.omics_projection = nn.Identity()
            self.esm2_projection = nn.Sequential(
                nn.Linear(self.esm2_dim, self.hidden_dim),
                nn.LayerNorm(self.hidden_dim),
                nn.GELU(),
                nn.Dropout(self.dropout),
            )
            self.gate_linear = nn.Linear(self.omics_dim + self.hidden_dim, self.hidden_dim)
            self.output_dim = self.omics_dim + self.hidden_dim

    def forward(self, x_omics: torch.Tensor, x_esm2: torch.Tensor):
        if self.fusion_mode == "gated":
            h_omics = self.omics_projection(x_omics)
            h_esm2 = self.esm2_projection(x_esm2)
            gate = torch.sigmoid(self.gate_linear(torch.cat([h_omics, h_esm2], dim=-1)))
            fused = gate * h_omics + (1.0 - gate) * h_esm2
            return fused, gate

        h_esm2 = self.esm2_projection(x_esm2)
        gate = torch.sigmoid(self.gate_linear(torch.cat([x_omics, h_esm2], dim=-1)))
        gated_esm2 = gate * h_esm2
        fused = torch.cat([x_omics, gated_esm2], dim=-1)
        return fused, gate
