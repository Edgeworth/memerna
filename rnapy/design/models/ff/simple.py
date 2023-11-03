from typing import Any

import torch
import torch.nn.functional as F
from torch import nn

from rnapy.design.harness.model import Model
from rnapy.design.rna.pipeline_cfg import RnaPipelineCfg


class SimpleFF(Model):
    """Simple feed forward network, e.g. to do a basic check the pipeline is working."""

    model: nn.Module
    cfg: RnaPipelineCfg

    def __init__(self, *, cfg: RnaPipelineCfg) -> None:
        super().__init__()
        self.cfg = cfg

        self.model = nn.Sequential(
            nn.Flatten(1),
            nn.Linear(cfg.max_seq_len * cfg.tensor.db_dim(), self.cfg.dim_feedforward),
            nn.ReLU(),
            nn.Linear(self.cfg.dim_feedforward, self.cfg.dim_feedforward),
            nn.ReLU(),
            nn.Linear(self.cfg.dim_feedforward, cfg.max_seq_len * cfg.tensor.primary_dim()),
            # Output (batch_size, d_primary, seq_len) so it fits nicely into cross entropy.
            nn.Unflatten(1, (cfg.tensor.primary_dim(), cfg.max_seq_len)),
        )

    def forward(self, db: torch.Tensor, _: torch.Tensor) -> Any:
        """
        Args:
            db: db structure, shape (batch_size, seq_len)

        Returns:
            output tensor, shape (batch_size, seq_len, d_primary)
        """

        # Convert to one hot encoding: (batch_size, seq_len, db_dim)
        db = F.one_hot(db, num_classes=self.d_emb).float()

        # Input tensor, shape (batch_size, seq_len, db_dim)
        return [self.model(db)]

    def model_prediction(self, out: torch.Tensor) -> torch.Tensor:
        return out.argmax(dim=-2)

    def model_loss(
        self, *, batch: list[torch.Tensor], outs: list[torch.Tensor]
    ) -> tuple[torch.Tensor, torch.Tensor]:
        out = outs[0]  # shape: (batch_size, seq_len, d_out_tok)
        y = batch[-1]  # primary, shape: (batch_size, seq_len)

        loss = F.cross_entropy(out, y)
        pred = self.model_prediction(out)

        return loss, (pred == y).type(torch.float)
