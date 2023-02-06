from typing import Any

from rnapy.design.harness.model import Model
from rnapy.design.rna.pipeline_cfg import RnaPipelineCfg
import torch
from torch import nn
import torch.nn.functional as F


# TODO(1): Finish.
class MlmFF(Model):
    """Simple feed forward network, but does masked language model stuff."""

    model: nn.Module
    cfg: RnaPipelineCfg

    def __init__(self, *, cfg: RnaPipelineCfg):
        super().__init__()
        self.cfg = cfg
        self.d_emb = max(cfg.tensor.db_dim(), cfg.tensor.primary_dim())

        self.model = nn.Sequential(
            nn.Flatten(1),
            nn.Linear(2 * cfg.max_seq_len * self.d_emb, self.cfg.dim_feedforward),
            nn.ReLU(),
            nn.Linear(self.cfg.dim_feedforward, self.cfg.dim_feedforward),
            nn.ReLU(),
            nn.Linear(self.cfg.dim_feedforward, cfg.max_seq_len * cfg.tensor.primary_dim()),
            # Output (batch_size, d_primary, seq_len) so it fits nicely into cross entropy.
            nn.Unflatten(1, (cfg.tensor.primary_dim(), cfg.max_seq_len)),
        )

    def forward(
        self,
        db: torch.Tensor,
        primary: torch.Tensor,
    ) -> Any:
        """
        Args:
            db: db structure, shape (batch_size, seq_len)
            primary: in primary structure, shape (batch_size, seq_len)

        Returns:
            output tensor, shape (batch_size, seq_len, d_primary)
        """
        # TODO(1): do masking.

        # Convert to one hot encoding: (batch_size, seq_len, d_emb)
        # TODO(0): undo db here, use masking? Or MLMFeedForward
        primary = F.one_hot(db, num_classes=self.d_emb).float()
        db = F.one_hot(db, num_classes=self.d_emb).float()

        # Input tensor, shape (batch_size, 2 * seq_len, d_emb)
        x = torch.cat([db, primary], dim=-2)
        return [self.model(x)]

    def model_prediction(self, out: torch.Tensor) -> torch.Tensor:
        return out.argmax(dim=-2)

    def model_loss(
        self,
        *,
        batch: list[torch.Tensor],
        outs: list[torch.Tensor],
    ) -> tuple[torch.Tensor, torch.Tensor]:
        out = outs[0]  # shape: (batch_size, seq_len, d_out_tok)
        y = batch[-1]  # primary, shape: (batch_size, seq_len)

        loss = F.cross_entropy(out, y)
        pred = self.model_prediction(out)

        return loss, (pred == y).type(torch.float)
