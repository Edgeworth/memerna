from typing import Any

from rnapy.design.harness.model import Model
from rnapy.design.rna.tensor import RnaTensor
from rnapy.design.transformer.transformer_model import TransformerModel
import torch
from torch import nn


class RnaTransformer(Model):
    model: TransformerModel

    def __init__(self, *, d_emb: int, max_seq_len: int):
        super().__init__()

        # Input here is the db structure, and the output is the primary structure.
        self.model = TransformerModel(
            d_seq=max_seq_len,
            d_inp_tok=RnaTensor.db_dim(),
            d_out_tok=RnaTensor.primary_dim(),
            d_emb=d_emb,
            dropout=0.1,
        )

    def forward(
        self,
        db: torch.Tensor,
        primary: torch.Tensor,
    ) -> Any:
        """
        Args:
            db: db structure, shape [batch_size, seq_len, d_db]
            primary: primary structure, shape [batch_size, seq_len, d_primary]

        Returns:
            output tensor, shape [batch_size, seq_len, d_out_tok]
        """

        device = db.device
        out_seq = primary[:, :-1]
        inp_mask = torch.zeros(db.shape[1], db.shape[1]).to(device)
        out_mask = nn.Transformer.generate_square_subsequent_mask(out_seq.shape[1]).to(device)

        return self.model(inp_seq=db, out_seq=out_seq, inp_mask=inp_mask, out_mask=out_mask)

    def model_prediction(self, out: torch.Tensor) -> torch.Tensor:
        return out.argmax(dim=-1)

    def model_loss(
        self,
        *,
        batch: list[torch.Tensor],
        out: torch.Tensor,
        loss_fn: nn.Module,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        y = batch[1][:, 1:]
        loss = loss_fn(out.reshape(-1, RnaTensor.primary_dim()), y.reshape(-1))
        correct = (self.model_prediction(out) == y).type(torch.float)
        return loss, correct
