from typing import Any

from rnapy.design.harness.model import Model
from rnapy.design.transformer.transformer_model import TransformerModel
import torch
from torch import nn

IN_TOKEN = 3
OUT_TOKEN = 5
EMB_SIZE = 128 # This can't be too small - need to be able to encode position.

class RnaTransformer(Model):
    model: TransformerModel

    def __init__(self, *, max_seq_len: int):
        super().__init__()

        # Input here is the secondary structure, and the output is the primary structure.
        self.model = TransformerModel(
            d_seq=max_seq_len,
            d_inp_tok=IN_TOKEN,  # TODO: abstract this based on mapping choice.
            d_out_tok=OUT_TOKEN,  # TODO: same here
            d_emb=EMB_SIZE,  # TODO: Parameter to adjust.
            dropout=0.1,
        )

    def forward(
        self,
        secondary: torch.Tensor,
        primary: torch.Tensor,
    ) -> Any:
        """
        Args:
            secondary: secondary structure, shape [batch_size, seq_len, d_secondary]
            primary: primary structure, shape [batch_size, seq_len, d_primary]

        Returns:
            output tensor, shape [batch_size, seq_len, d_out_tok]
        """

        device = secondary.device
        inp_mask = torch.zeros(secondary.shape[1], secondary.shape[1]).to(device)
        out_mask = nn.Transformer.generate_square_subsequent_mask(primary.shape[1]).to(device)

        return self.model(inp_seq=secondary, out_seq=primary, inp_mask=inp_mask, out_mask=out_mask)

    def model_prediction(self, out: torch.Tensor) -> torch.Tensor:
        return out.argmax(dim=-1)

    def model_loss(
        self,
        *,
        batch: list[torch.Tensor],
        out: torch.Tensor,
        loss_fn: nn.Module,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        y = batch[1]
        loss = loss_fn(out.reshape(-1, OUT_TOKEN), y.reshape(-1))
        correct = (self.model_prediction(out) == y).type(torch.float)
        return loss, correct
