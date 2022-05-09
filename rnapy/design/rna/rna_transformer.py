from math import ceil
from typing import Any

from rnapy.design.harness.model import Model
from rnapy.design.rna.config import RnaPipelineConfig
from rnapy.design.rna.tensor import MASK_IDX, PAD_IDX, RnaTensor
from rnapy.design.transformer.transformer_model import TransformerModel
import torch
import torch.nn.functional as F
import random


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
    def mask_for_values(src: torch.Tensor, values: list[int]) -> torch.Tensor:
        """Returns a mask tensor of the same shape as src, true iff it is one of the given values"""
        mask = torch.full_like(src, False, dtype=torch.bool)
        for v in values:
            mask |= src == v
        return mask.to(src.device)

    @staticmethod
    def prob_mask_subset(mask: torch.Tensor, prob: float) -> torch.Tensor:
        """Returns mask that is a subset of the given mask. The last dimension will have a
        subset selected according to prob. This avoids some rows being entirely empty or full."""
        # Get each row of the mask (last dimension)
        new_mask = mask.detach()
        new_mask = new_mask.reshape((-1, new_mask.shape[-1]))
        mask_list = new_mask.tolist()
        for row in mask_list:
            cur_num_selected = sum(int(i) for i in row)
            num_subset = int(ceil(cur_num_selected * prob))
            index_map = []
            cur_idx = 0
            for i, v in enumerate(row):
                if v:
                    index_map.append(i)
                cur_idx += 1

            # Reset all row values to False
            for i in range(len(row)):
                row[i] = False

            # Grab indices into the currently selected masked values without replacement.
            for index in random.sample(range(cur_num_selected), num_subset):
                row[index_map[index]] = True
        return torch.BoolTensor(mask_list).reshape(mask.shape).to(mask.device)

    @staticmethod
    def prob_mask_like(src: torch.Tensor, prob: float) -> torch.Tensor:
        return RnaTransformer.prob_mask_subset(torch.ones_like(src, dtype=torch.bool), prob)

    def forward(
        self,
        db: torch.Tensor,
        primary: torch.Tensor,
        randomize_mask: bool = True,
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
        # If not randomizing the mask ourselves, input should come with masked values already.
        mask = torch.zeros_like(primary, dtype=torch.bool)
        if randomize_mask:
            # Decide which bases to mask. First find all non-padded tokens.
            mask = ~self.mask_for_values(primary, [PAD_IDX])

            # Select subset with given probability to mask out.
            mask = self.prob_mask_subset(mask, self.cfg.mask_prop)

            # TODO: Add purposely incorrect ones.
            # Select some of the masked tokens to be passed through as correct.
            passthrough_mask = self.prob_mask_like(mask, self.cfg.mask_passthrough_prop)

            # Tokens to be actually masked out
            blocked_mask = mask & ~passthrough_mask

            out_seq = out_seq.masked_fill(blocked_mask, MASK_IDX)

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
