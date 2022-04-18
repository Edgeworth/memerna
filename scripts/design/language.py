from pathlib import Path
from typing import Any, Iterator

import click
from scripts.design.trainer import Trainer
from torchtext.data.utils import get_tokenizer
from torchtext.vocab import build_vocab_from_iterator, Vocab
from torch.utils.data import Dataset, ConcatDataset
import torch
import torchtext

from scripts.design.transformer.transformer_model import TransformerModel

MAX_SEQ_LEN = 32  # max sequence length


class SequenceDataset(Dataset):
    """Dataset that produces overlapped sequences from the given data and predicts the next value"""

    seq: torch.Tensor
    seq_size: int

    def __init__(self, *, seq: torch.Tensor, seq_size: int) -> None:
        if not isinstance(seq, torch.Tensor):
            raise TypeError("seq must be a torch.Tensor")
        self.seq = seq
        self.seq_size = seq_size

    def __len__(self) -> int:
        return max(len(self.seq) - self.seq_size, 0)

    def __getitem__(self, idx: int) -> tuple[torch.Tensor, Any]:
        """
        Returns:
            Next value in the sequence for each point in the sequence"""
        return (self.seq[idx : idx + self.seq_size], self.seq[idx + 1 : idx + self.seq_size + 1])


class LanguagePipeline:
    output_path: Path
    train_data: Dataset
    valid_data: Dataset
    model: TransformerModel
    trainer: Trainer
    vocab: Vocab

    def __init__(self, *, output_path: Path) -> None:
        self.output_path = output_path

        train_data = torchtext.datasets.WikiText2(
            root=output_path,
            split="train",
        )
        self._build_vocab(train_data)

        self.train_data = self._build_data(train_data)

        self.valid_data = self._build_data(
            torchtext.datasets.WikiText2(
                root=output_path,
                split="valid",
            )
        )

        self.model = TransformerModel(
            d_seq=MAX_SEQ_LEN,
            d_inp_words=len(self.vocab),
            d_out_words=len(self.vocab),
            d_emb=256,
            dropout=0.1,
        )
        self.trainer = Trainer(
            model=self.model,
            train_data=self.train_data,
            valid_data=self.valid_data,
            path=output_path,
            clip_grad_norm=1.0,
        )

    def _build_data(self, data: Iterator[str]) -> Dataset:
        tokenizer = get_tokenizer("basic_english")
        datasets = []
        for text in data:
            tokens = torch.LongTensor(self.vocab(tokenizer(text)))
            if len(tokens) <= MAX_SEQ_LEN:
                continue
            datasets.append(SequenceDataset(seq=tokens, seq_size=MAX_SEQ_LEN))
        return ConcatDataset(datasets)

    def _build_vocab(self, data: Iterator[str]) -> None:
        tokenizer = get_tokenizer("basic_english")
        self.vocab = build_vocab_from_iterator(map(tokenizer, data), specials=["<unk>"])
        self.vocab.set_default_index(self.vocab["<unk>"])

    def run(self, epochs: int) -> None:
        click.echo(f"Running transformer at {self.output_path}")
        self.trainer.run(epochs)
