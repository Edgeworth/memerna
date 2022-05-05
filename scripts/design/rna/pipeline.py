from pathlib import Path

import click
from scripts.design.harness.config import TrainConfig
from scripts.design.harness.trainer import Trainer
from scripts.design.rna.dataset import RnaDataset
from scripts.design.transformer.transformer_model import TransformerModel
from torch.utils.data import Dataset

MAX_SEQ_LEN = 32  # max sequence length


class RnaPipeline:
    output_path: Path
    train_data: Dataset
    valid_data: Dataset
    test_data: Dataset
    trainer: Trainer

    def __init__(self, *, output_path: Path, checkpoint_path: Path | None) -> None:
        self.output_path = output_path

        self.train_data = RnaDataset(num_seq=1000, seq_size=128)
        self.valid_data = RnaDataset(num_seq=1000, seq_size=128)
        self.test_data = RnaDataset(num_seq=1000, seq_size=128)

        model = TransformerModel(
            d_seq=MAX_SEQ_LEN,
            d_inp_word=5,  # TODO: make class for rna one hot vec
            d_out_word=5,
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

    # def _build_data(self, data: Iterator[str]) -> Dataset:
    #     tokenizer = get_tokenizer("basic_english")
    #     datasets = []
    #     for text in data:
    #         tokens = torch.LongTensor(self.vocab(tokenizer(text)))
    #         if len(tokens) <= MAX_SEQ_LEN:
    #             continue
    #         datasets.append(SequenceDataset(seq=tokens, seq_size=MAX_SEQ_LEN))
    #     return ConcatDataset(datasets)

    # def _indices_to_words(self, indices: torch.Tensor) -> list[str]:
    #     return self.vocab.lookup_tokens(indices.tolist())

    # def _words_to_indices(self, words: list[str]) -> torch.Tensor:
    #     return torch.LongTensor(self.vocab(words))

    # def _predict_next(
    #     self,
    #     words: list[str],
    # ) -> list[str]:
    #     X = self._words_to_indices(words)[-MAX_SEQ_LEN:]
    #     X = X.unsqueeze(0)
    #     dm = self.trainer.optimizer.dm
    #     out = dm(X)
    #     pred = dm.prediction(out=out).to("cpu")
    #     pred_words = self._indices_to_words(pred.squeeze(0))
    #     return words + [pred_words[-1]]

    def predict(self, start: str) -> None:
        # Choose random element from the test set:
        print("Start input: ", start)
        # for _ in range(50):
        #     words = self._predict_next(words)
        # print(" ".join(words))

    def train(self, epochs: int) -> None:
        click.echo(f"Running transformer at {self.output_path}")
        self.trainer.run(epochs)
