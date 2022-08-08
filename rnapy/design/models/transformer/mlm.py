from math import ceil
import random
from typing import Any

from rnapy.design.harness.model import Model
from rnapy.design.models.parts.transformer.transformer_model import TransformerModel
from rnapy.design.rna.pipeline_cfg import RnaPipelineCfg
from rnapy.design.rna.tensor import MASK_IDX
from rnapy.design.rna.tensor import PAD_IDX
import torch
import torch.nn.functional as F


class MLMTransformer(Model):
    model: TransformerModel
    cfg: RnaPipelineCfg

    def __init__(self, *, cfg: RnaPipelineCfg):
        super().__init__()
        self.cfg = cfg

        # Input here is the db structure, and the output is the primary structure.
        self.model = TransformerModel(
            d_seq=cfg.max_seq_len,
            d_inp_tok=cfg.tensor.db_dim(),
            d_out_tok=cfg.tensor.primary_dim(),
            d_emb=cfg.mlm.d_emb,
            dropout=0.1,
            nhead=cfg.mlm.nhead,
            num_encoder_layers=cfg.mlm.num_encoder_layers,
            num_decoder_layers=cfg.mlm.num_decoder_layers,
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
        return MLMTransformer.prob_mask_subset(torch.ones_like(src, dtype=torch.bool), prob)

    def randomize_mask(self, primary: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        # Tokens to not randomize.
        ignored_tokens = [PAD_IDX]

        # Decide which bases to mask. First find all non-padded tokens.
        mask = ~self.mask_for_values(primary, ignored_tokens)

        # Select subset with given probability to mask out.
        mask = self.prob_mask_subset(mask, self.cfg.mlm.mask_prop)

        # TODO(1): Add purposely incorrect ones.
        # Select some of the masked tokens to be passed through as correct.
        passthrough_prop = self.cfg.mlm.mask_correct_prop + self.cfg.mlm.mask_incorrect_prop
        passthrough_mask = self.prob_mask_subset(mask, passthrough_prop)

        # Tokens to be actually masked out
        blocked_mask = mask & ~passthrough_mask

        primary = primary.masked_fill(blocked_mask, MASK_IDX)

        if self.cfg.mlm.mask_incorrect_prop > 0:
            # Select some of the masked tokens to be passed through as incorrect.
            random_tokens = torch.randint(
                0,
                self.cfg.tensor.primary_dim(),
                size=primary.shape,
                device=primary.device,
            )
            # Remove the ignored tokens. We may not passthrough enough incorrect tokens
            # but oh well.
            random_ignored = self.mask_for_values(random_tokens, ignored_tokens)
            incorrect_mask = self.prob_mask_subset(
                passthrough_mask,
                self.cfg.mlm.mask_incorrect_prop,
            )
            incorrect_mask &= ~random_ignored
            primary = torch.where(incorrect_mask, random_tokens, primary)

        return (primary, mask)

    def forward(
        self,
        db: torch.Tensor,
        primary: torch.Tensor,
        randomize_mask: bool = True,
    ) -> Any:
        """
        Args:
            db: db structure, shape (batch_size, seq_len, d_db)
            primary: in primary structure, shape (batch_size, seq_len, d_primary)
            randomize_mask: Whether to make a random mask (used during training).

        Returns:
            output tensor, shape (batch_size, seq_len, d_out_tok)
            mask tensor (if randomly generated), shape (batch_size, seq_len, d_primary)
        """

        # If not randomizing the mask ourselves, input should come with masked values already.
        mask = torch.zeros_like(primary, dtype=torch.bool)
        if randomize_mask:
            primary, mask = self.randomize_mask(primary)

        # print(primary[0])
        # print(self.cfg.tensor.to_primary(primary[0]))
        # print(self.cfg.tensor.to_db(db[0]))
        return [self.model(inp_seq=db, out_seq=primary), mask]

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

        # print(self.cfg.tensor.to_primary(y[0]))
        # print(self.cfg.tensor.to_primary(self.model_prediction(out)[0]))
        return loss, accuracy
