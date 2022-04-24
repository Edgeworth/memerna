from typing import Any

from scripts.design.harness.model import Model
from scripts.design.transformer.positional_encoder import PositionalEncoder
from scripts.design.transformer.word_embedding import WordEmbedding
import torch
from torch import nn


class TransformerModel(Model):
    d_out_words: int

    def __init__(
        self,
        *,
        d_seq: int,
        d_inp_word: int,
        d_out_word: int,
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
        self.d_out_words = d_out_word

        self.inp_emb = WordEmbedding(d_word=d_inp_word, d_emb=d_emb)
        self.out_emb = WordEmbedding(d_word=d_out_word, d_emb=d_emb)
        self.pos_encoder = PositionalEncoder(d_emb=d_emb, max_seq_len=d_seq, dropout=dropout)
        # TODO: Think about parameters here.
        self.transformer = nn.Transformer(
            d_model=d_emb,
            nhead=8,  # TODO: think about this.
            num_encoder_layers=6,  # TODO: think about this.
            num_decoder_layers=6,  # TODO: think about this.
            dim_feedforward=64,  # TODO: think about this.
            batch_first=True,
            dropout=dropout,
        )
        # Take output of transformer and pass through linear layer to predict a
        # word in the output vocabulary.
        self.linear = nn.Linear(d_emb, d_out_word)

    def forward(
        self,
        inp_seq: torch.Tensor,
        out_seq: torch.Tensor,
        inp_mask: torch.Tensor,
        out_mask: torch.Tensor,
    ) -> Any:
        """
        Args:
            inp: input tensor of sequence of "words" by their indices, shape [batch_size, seq_len]
            out: output tensor of sequence of "words" by their indices, shape [batch_size, seq_len]
            inp_mask: mask of input sequence, shape [seq_len, seq_len]
            out_mask: mask of output sequence, shape [seq_len, seq_len]

        Returns:
            output tensor, shape [batch_size, seq_len, d_out_word]
        """
        # Embedding for input, output shape: [batch_size, seq_len, d_emb]
        inp_seq = self.pos_encoder(self.inp_emb(inp_seq))
        # Embedding for output, ouptut shape: [batch_size, seq_len, d_emb]
        out_seq = self.pos_encoder(self.out_emb(out_seq))

        # Transformer, output shape: [batch_size, seq_len, d_emb]
        attn = self.transformer(src=inp_seq, tgt=out_seq, src_mask=inp_mask, tgt_mask=out_mask)

        # Linear layer, output shape: [batch_size, seq_len, d_out_word]
        out = self.linear(attn)

        return out

    def model_inputs(self, X: torch.Tensor) -> list[Any]:
        # Input mask adds 0's, so no effect and no masking.
        inp_mask = torch.zeros(X.shape[1], X.shape[1])
        # Output mask adds an upper triangular matrix of -inf, so
        # the softmax'd outputs for them are zero.
        out_mask = nn.Transformer.generate_square_subsequent_mask(X.shape[1])
        return [X, X, inp_mask, out_mask]

    def model_prediction(self, out: torch.Tensor) -> torch.Tensor:
        return out.argmax(dim=-1)

    def model_loss(
        self,
        *,
        out: torch.Tensor,
        y: torch.Tensor,
        loss_fn: nn.Module,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        loss = loss_fn(out.reshape(-1, self.d_out_words), y.reshape(-1))
        correct = (self.model_prediction(out) == y).type(torch.float)
        return loss, correct
