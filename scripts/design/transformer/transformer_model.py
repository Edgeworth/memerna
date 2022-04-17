from typing import Any

from scripts.design.transformer.positional_encoder import PositionalEncoder
from scripts.design.transformer.word_embedding import WordEmbedding
import torch
from torch import nn


class TransformerModel(nn.Module):
    def __init__(
        self,
        *,
        d_seq: int,
        d_inp_words: int,
        d_out_words: int,
        d_emb: int,
        dropout: float = 0.1,
    ) -> None:
        """
        Args:
            d_seq: dimension of input/output sequence (input/output sequence length)
            d_inp_words: dimension of input words
            d_out_words: dimension of output words
            d_emb: dimension of embedding of words
            dropout: dropout rate
        """
        super().__init__()
        self.inp_emb = WordEmbedding(d_words=d_inp_words, d_emb=d_emb)
        self.out_emb = WordEmbedding(d_words=d_out_words, d_emb=d_emb)
        self.pos_encoder = PositionalEncoder(d_emb=d_emb, max_seq_len=d_seq, dropout=dropout)
        self.transformer = nn.Transformer(
            d_model=d_emb,
            nhead=8,
            num_encoder_layers=6,
            num_decoder_layers=6,
            dim_feedforward=2048,
            batch_first=True,
            dropout=dropout,
        )
        # Take output of transformer and pass through linear layer to predict a
        # word in the output vocabulary.
        self.linear = nn.Linear(d_emb, d_out_words)

    def forward(self, inp_seq: torch.Tensor, out_seq: torch.Tensor) -> Any:
        """
        Args:
            inp: input tensor of sequence of "words", shape [batch_size, seq_len, d_word]
            out: output tensor of sequence of "words", shape [batch_size, seq_len, d_word]

        Returns:
            , shape [batch_size, seq_len, d_out_word]
        """
        # Embedding for input, shape: [batch_size, seq_len, d_emb]
        inp_seq = self.pos_encoder(self.inp_emb(inp_seq))
        # Embedding for output, shape: [batch_size, seq_len, d_emb]
        out_seq = self.pos_encoder(self.out_emb(out_seq))

        # Transformer, shape: [batch_size, seq_len, d_emb]
        out = self.transformer(src=inp_seq, tgt=out_seq)
        return self.linear(out)
