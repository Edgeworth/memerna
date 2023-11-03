from typing import Any

import torch
from torch import nn

from rnapy.design.models.parts.transformer.positional_encoder import PositionalEncoder
from rnapy.design.models.parts.transformer.token_embedding import TokenEmbedding


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
        nhead: int = 8,
        num_encoder_layers: int = 6,
        num_decoder_layers: int = 6,
        dim_feedforward: int = 2048,
        activation: str = "relu",
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
        self.transformer = nn.Transformer(
            batch_first=True,
            d_model=d_emb,
            nhead=nhead,
            num_encoder_layers=num_encoder_layers,
            num_decoder_layers=num_decoder_layers,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            activation=activation,
        )
        # Take output of transformer and pass through linear layer to predict a
        # token in the output vocabulary.
        self.linear = nn.Linear(d_emb, d_out_tok)

    def forward(
        self,
        inp_seq: torch.Tensor,
        out_seq: torch.Tensor,
        inp_mask: torch.Tensor | None = None,
        out_mask: torch.Tensor | None = None,
        inp_token_mask: torch.ByteTensor | torch.BoolTensor | None = None,
        out_token_mask: torch.ByteTensor | torch.BoolTensor | None = None,
    ) -> Any:
        """
        Args:
            inp: input tensor of sequence of "tokens" by their indices, shape (batch_size, seq_len)
            out: output tensor of sequence of "tokens" by their indices, shape (batch_size, seq_len)
            inp_mask: mask of input attention weights, shape (seq_len, seq_len)
            out_mask: mask of output attention weights, shape (seq_len, seq_len)
            inp_token_mask: non-zero to ignore tokens in input, shape ([batch_size], seq_len)
            out_token_mask: non-zero to ignore tokens in output, shape ([batch_size], seq_len)

        Returns:
            output tensor, shape (batch_size, seq_len, d_out_tok)
        """

        # Embedding for input, output shape: [batch_size, seq_len, d_emb]
        inp_seq = self.pos_encoder(self.inp_emb(inp_seq))
        # Embedding for output, ouptut shape: [batch_size, seq_len, d_emb]
        out_seq = self.pos_encoder(self.out_emb(out_seq))

        # Transformer, output shape: [batch_size, seq_len, d_emb]
        attn = self.transformer(
            src=inp_seq,
            tgt=out_seq,
            src_mask=inp_mask,
            tgt_mask=out_mask,
            src_key_padding_mask=inp_token_mask,
            tgt_key_padding_mask=out_token_mask,
        )

        # Linear layer, output shape: [batch_size, seq_len, d_out_tok]
        return self.linear(attn)
