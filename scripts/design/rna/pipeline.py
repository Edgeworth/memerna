from pathlib import Path
from typing import Iterator

import click
from scripts.design.harness.config import TrainConfig
from scripts.design.harness.sequence import SequenceDataset
from scripts.design.harness.trainer import Trainer
from scripts.design.transformer.transformer_model import TransformerModel
import torch
from torch.utils.data import ConcatDataset
from torch.utils.data import Dataset
import torchtext
from torchtext.data.utils import get_tokenizer
from torchtext.vocab import build_vocab_from_iterator
from torchtext.vocab import Vocab

MAX_SEQ_LEN = 32  # max sequence length


class RnaPipeline:
    output_path: Path
    train_data: Dataset
    valid_data: Dataset
    test_data: Dataset
    trainer: Trainer
    vocab: Vocab

    def __init__(self, *, output_path: Path, checkpoint_path: Path | None) -> None:
        self.output_path = output_path

        train_data = torchtext.datasets.WikiText2(
            root=output_path,
            split="train",
        )
        self._build_vocab(train_data)

        self.train_data = self._build_data(train_data)
        self.valid_data = self._build_data(
            torchtext.datasets.WikiText2(root=output_path, split="valid"),
        )
        self.test_data = self._build_data(
            torchtext.datasets.WikiText2(root=output_path, split="test"),
        )

        model = TransformerModel(
            d_seq=MAX_SEQ_LEN,
            d_inp_word=len(self.vocab),
            d_out_word=len(self.vocab),
            d_emb=512,  # TODO: Parameter to adjust.
            dropout=0.1,
        )
        cfg = TrainConfig(
            model_name="LanguageTransformer",
            output_path=output_path,
            profile=False,
            batch_size=32,
            train_samples=10000,
            fast_valid_samples=512,
            accurate_valid_samples=4096,
            clip_grad_norm=1.0,
        )
        self.trainer = Trainer(
            model=model,
            train_data=self.train_data,
            valid_data=self.valid_data,
            cfg=cfg,
            checkpoint_path=checkpoint_path,
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

    def _indices_to_words(self, indices: torch.Tensor) -> list[str]:
        return self.vocab.lookup_tokens(indices.tolist())

    def _words_to_indices(self, words: list[str]) -> torch.Tensor:
        return torch.LongTensor(self.vocab(words))

    def _predict_next(
        self,
        words: list[str],
    ) -> list[str]:
        X = self._words_to_indices(words)[-MAX_SEQ_LEN:]
        X = X.unsqueeze(0)
        dm = self.trainer.optimizer.dm
        out = dm(X)
        pred = dm.prediction(out=out).to("cpu")
        pred_words = self._indices_to_words(pred.squeeze(0))
        return words + [pred_words[-1]]

    def predict(self, start: str) -> None:
        # Choose random element from the test set:
        print("Start input: ", start)
        tokenizer = get_tokenizer("basic_english")
        words = tokenizer(start)
        for _ in range(50):
            words = self._predict_next(words)
        print(" ".join(words))

    def train(self, epochs: int) -> None:
        click.echo(f"Running transformer at {self.output_path}")
        self.trainer.run(epochs)
