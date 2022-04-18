from pathlib import Path

import click
from scripts.design.trainer import Trainer
from torchtext.data.utils import get_tokenizer
from torchtext.vocab import build_vocab_from_iterator, Vocab
import torchtext

from scripts.design.transformer.transformer_model import TransformerModel

SEQ_LEN =


class LanguagePipeline:
    output_path: Path
    train_data: torchtext.datasets.WikiText2
    valid_data: torchtext.datasets.WikiText2
    model: TransformerModel
    trainer: Trainer
    vocab: Vocab

    def __init__(self, *, output_path: Path) -> None:
        self.output_path = output_path

        self.train_data = torchtext.datasets.WikiText2(
            root=output_path,
            split="train",
        )
        self.valid_data = torchtext.datasets.WikiText2(
            root=output_path,
            split="valid",
        )

        print(self.train_data.iter().next())

        self._build_vocab()

        self.model = TransformerModel(
        )
        self.trainer = Trainer(
            model=self.model,
            train_data=self.train_data,
            valid_data=self.valid_data,
            path=output_path,
            clip_grad_norm=1.0,
        )

    def _build_vocab(self) -> None:
        tokenizer = get_tokenizer("basic_english")
        vocab = build_vocab_from_iterator(map(tokenizer, self.train_data), specials=["<unk>"])
        vocab.set_default_index(vocab["<unk>"])

    def run(self, epochs: int) -> None:
        click.echo(f"Running transformer at {self.output_path}")

        self.trainer.run(epochs)
