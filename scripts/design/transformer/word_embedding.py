import math
from typing import Any

import torch
from torch import nn


class WordEmbedding(nn.Module):
    def __init__(self, *, d_words: int, d_emb: int) -> None:
        """
        Args:
            d_words: dimension of vocabulary (size of vocabulary)
            d_emb: dimension of embedding of words
        """
        super().__init__()
        self.emb = nn.Embedding(d_words, d_emb)
        self.d_emb = d_emb

    def forward(self, seq: torch.Tensor) -> Any:
        # Transformers multiply the embedding vector by sqrt(d_emb).
        return self.emb(seq.long()) * math.sqrt(self.d_emb)
