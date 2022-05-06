from typing import Any

from rnapy.design.transformer.positional_encoder import PositionalEncoder
from rnapy.design.transformer.token_embedding import TokenEmbedding
import torch
from torch import nn


class TransformerModel(nn.Module):
    """Transformer model. Good to composte this into a Model subclass."""

    def __init__(
        self,
        *,
        d_seq: int,
        d_inp_tok: int,
        d_out_tok: int,
        d_emb: int,
        dropout: float = 0.1,
    ) -> None:
        """
        Args:
            d_seq: dimension of input/output sequence (input/output sequence length)
            d_inp_tok: dimension of input tokens
            d_out_tok: dimension of output tokens
            d_emb: dimension of embedding of tokens
            dropout: dropout rate
        """
        super().__init__()

        self.inp_emb = TokenEmbedding(d_tok=d_inp_tok, d_emb=d_emb)
        self.out_emb = TokenEmbedding(d_tok=d_out_tok, d_emb=d_emb)
        self.pos_encoder = PositionalEncoder(d_emb=d_emb, max_seq_len=d_seq, dropout=dropout)
        # TODO: Think about parameters here.
        self.transformer = nn.Transformer(
            d_model=d_emb,
            nhead=8,  # TODO: Parameter to adjust.
            num_encoder_layers=6,  # TODO: Parameter to adjust.
            num_decoder_layers=6,  # TODO: Parameter to adjust.
            dim_feedforward=2048,  # TODO: Parameter to adjust.
            batch_first=True,
            dropout=dropout,
        )
        # Take output of transformer and pass through linear layer to predict a
        # token in the output vocabulary.
        self.linear = nn.Linear(d_emb, d_out_tok)

    def forward(
        self,
        inp_seq: torch.Tensor,
        out_seq: torch.Tensor,
        inp_mask: torch.Tensor,
        out_mask: torch.Tensor,
    ) -> Any:
        """
        Args:
            inp: input tensor of sequence of "tokens" by their indices, shape [batch_size, seq_len]
            out: output tensor of sequence of "tokens" by their indices, shape [batch_size, seq_len]
            inp_mask: mask of input sequence, shape [seq_len, seq_len]
            out_mask: mask of output sequence, shape [seq_len, seq_len]

        Returns:
            output tensor, shape [batch_size, seq_len, d_out_tok]
        """

        # Embedding for input, output shape: [batch_size, seq_len, d_emb]
        inp_seq = self.pos_encoder(self.inp_emb(inp_seq))
        # Embedding for output, ouptut shape: [batch_size, seq_len, d_emb]
        out_seq = self.pos_encoder(self.out_emb(out_seq))

        # Transformer, output shape: [batch_size, seq_len, d_emb]
        attn = self.transformer(src=inp_seq, tgt=out_seq, src_mask=inp_mask, tgt_mask=out_mask)

        # Linear layer, output shape: [batch_size, seq_len, d_out_tok]
        out = self.linear(attn)

        return out
