import math
from typing import Any

import torch
from torch import nn


class PositionalEncoder(nn.Module):
    """Positional encoding for transformer."""

    pos_emb: torch.Tensor

    def __init__(self, *, d_emb: int, max_seq_len: int, dropout: float = 0.1) -> None:
        """

        Args:
            d_emb: dimension of embedding of input
            max_seq_len: max length of sequence of input
            dropout: dropout rate
        """
        super().__init__()

        # Compute 1 / 10000 ^ (2i / d_model) to use as input for sin and cos.
        div = torch.exp(-torch.arange(0, d_emb, 2) * math.log(10000) / d_emb)

        # Compute position index (index of value in input sequence).
        pos = torch.arange(max_seq_len).unsqueeze(1)

        # Compute the positional encoding to add on to the input embedding.
        pos_emb = torch.zeros(max_seq_len, d_emb)
        pos_emb[:, 0::2] = torch.sin(pos * div)
        pos_emb[:, 1::2] = torch.cos(pos * div)

        self.dropout = nn.Dropout(dropout)

        # Positional embedding isn't a model parameter, but we want to save it
        # in the model, so register it as a buffer.
        self.register_buffer("pos_emb", pos_emb)

    def forward(self, seq_emb: torch.Tensor) -> Any:
        """
        Args:
            seq_emb: embedding of input sequence, shape [batch_size, seq_len, d_emb]
        Returns:
            position encoded tensor, shape [batch_size, seq_len, d_emb]
        """
        return self.dropout(seq_emb + self.pos_emb[: seq_emb.size(1)])
