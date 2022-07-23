import math
from typing import Any

import torch
from torch import nn


class TokenEmbedding(nn.Module):
    def __init__(self, *, d_tok: int, d_emb: int) -> None:
        """
        Args:
            d_tok: dimension of vocabulary (size of vocabulary)
            d_emb: dimension of embedding of tokens
        """
        super().__init__()
        self.emb = nn.Embedding(d_tok, d_emb)
        self.d_emb = d_emb

    def forward(self, seq: torch.Tensor) -> Any:
        # Transformers multiply the embedding vector by sqrt(d_emb) for
        # some unclear reason.
        return self.emb(seq.long()) * math.sqrt(self.d_emb)
