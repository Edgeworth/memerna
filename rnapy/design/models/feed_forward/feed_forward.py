from typing import Any

from rnapy.design.harness.model import Model
from rnapy.design.rna.config import RnaPipelineConfig
from rnapy.design.rna.tensor import PAD_IDX
import torch
import torch.nn as nn


class SimpleFeedForward(Model):
    """Simple feed forward network, e.g. to do a basic check the pipeline is working."""

    model: nn.Module
    cfg: RnaPipelineConfig

    def __init__(self, *, cfg: RnaPipelineConfig):
        super().__init__()
        self.cfg = cfg

        self.model = nn.Sequential(
            nn.Linear(cfg.tensor.primary_dim() + cfg.tensor.db_dim(), self.cfg.dim_feedforward),
            nn.ReLU(),
            nn.Linear(self.cfg.dim_feedforward, self.cfg.dim_feedforward),
            nn.ReLU(),
            nn.Linear(self.cfg.dim_feedforward, cfg.tensor.primary_dim()),
        )

    def forward(
        self,
        db: torch.Tensor,
        primary: torch.Tensor,
    ) -> Any:
        """
        Args:
            db: db structure, shape (batch_size, seq_len, d_db)
            primary: in primary structure, shape (batch_size, seq_len, d_primary)

        Returns:
            output tensor, shape (batch_size, seq_len, d_primary)
        """

        # Input tensor, shape (batch_size, seq_len, d_db + d_primary)
        x = torch.cat([db, primary], dim=-1)
        return [self.model(x)]

    def model_prediction(self, out: torch.Tensor) -> torch.Tensor:
        return out.argmax(dim=-1)

    def model_loss(
        self,
        *,
        batch: list[torch.Tensor],
        outs: list[torch.Tensor],
    ) -> tuple[torch.Tensor, torch.Tensor]:
        out = outs[0]  # shape: (batch_size, seq_len, d_out_tok)
        mask = outs[1]

        # Blank out all non-masked tokens (ones we aren't trying to predict)
        # so they aren't included in the loss.
        y = batch[-1].masked_fill(~mask, PAD_IDX)

        loss = F.cross_entropy(
            out.reshape(-1, self.cfg.tensor.primary_dim()),
            y.reshape(-1),
            ignore_index=PAD_IDX,
        )

        # Select only the masked tokens for accuracy calculation.
        # This will collapse to 1D but that's okay.
        accuracy = (self.model_prediction(out) == y).type(torch.float).masked_select(mask)

        return loss, accuracy
