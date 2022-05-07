from typing import Any

from rnapy.design.harness.model import Model
from rnapy.design.rna.config import RnaPipelineConfig
from rnapy.design.rna.tensor import MASK_IDX, PAD_IDX, RnaTensor
from rnapy.design.transformer.transformer_model import TransformerModel
import torch
import torch.nn.functional as F


class RnaTransformer(Model):
    model: TransformerModel
    cfg: RnaPipelineConfig

    def __init__(self, *, cfg: RnaPipelineConfig):
        super().__init__()
        self.cfg = cfg

        # Input here is the db structure, and the output is the primary structure.
        self.model = TransformerModel(
            d_seq=cfg.max_seq_len,
            d_inp_tok=RnaTensor.db_dim(),
            d_out_tok=RnaTensor.primary_dim(),
            d_emb=cfg.d_emb,
            dropout=0.1,
            nhead=cfg.nhead,
            num_encoder_layers=cfg.num_encoder_layers,
            num_decoder_layers=cfg.num_decoder_layers,
            dim_feedforward=cfg.dim_feedforward,
            activation=cfg.activation,
        )

    @staticmethod
    def mask_values(src: torch.Tensor, values: list[int]) -> torch.Tensor:
        """Returns a mask tensor of the same shape as src, true iff it is one of the given values"""
        mask = torch.full_like(src, False, dtype=torch.bool)
        for v in values:
            mask |= src == v
        return mask

    @staticmethod
    def prob_mask_like(src: torch.Tensor, prob: float) -> torch.Tensor:
        return torch.rand_like(src, dtype=torch.float) < prob

    def forward(
        self,
        db: torch.Tensor,
        primary: torch.Tensor,
    ) -> Any:
        """
        Args:
            db: db structure, shape (batch_size, seq_len, d_db)
            primary_in: in primary structure, shape (batch_size, seq_len, d_primary)
            primary_out: expected primary structure, shape (batch_size, seq_len, d_primary)

        Returns:
            output tensor, shape (batch_size, seq_len, d_out_tok)
        """

        out_seq = primary
        # If not training, input should come with masked values already.
        mask = torch.zeros_like(primary, dtype=torch.bool)
        if self.training:
            # Decide which bases to mask. First find all non-padded tokens.
            mask = ~self.mask_values(primary, [PAD_IDX])

            # Select subset with given probability to mask out.
            mask = mask & self.prob_mask_like(mask, self.cfg.mask_prop)
            out_seq = out_seq.masked_fill(mask, MASK_IDX)

        return [self.model(inp_seq=db, out_seq=out_seq), mask]

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

        print(mask[0])

        # Blank out all non-masked tokens (ones we aren't trying to predict)
        # so they aren't included in the loss.
        y = batch[-1].masked_fill(~mask, PAD_IDX)

        loss = F.cross_entropy(
            out.reshape(-1, RnaTensor.primary_dim()), y.reshape(-1), ignore_index=PAD_IDX
        )

        # Select only the masked tokens for accuracy calculation.
        # Thils will collapse to 1D but that's okay.
        accuracy = (self.model_prediction(out) == y).type(torch.float).masked_select(mask)

        # print(RnaTensor.to_primary(batch[1][0]))
        # print(RnaTensor.to_primary(y[0]))
        # print(RnaTensor.to_primary(self.model_prediction(out)[0]))
        return loss, accuracy
