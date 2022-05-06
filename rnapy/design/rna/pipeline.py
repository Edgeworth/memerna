from pathlib import Path

import click
from rnapy.design.harness.config import TrainConfig
from rnapy.design.harness.trainer import Trainer
from rnapy.design.rna.dataset import RnaDataset
from rnapy.design.rna.rna_transformer import RnaTransformer
from torch.utils.data import Dataset
from rnapy.design.rna.tensor import BOS_IDX, RnaTensor
from rnapy.model.parse.sequence import db_to_secondary
import torch


class RnaPipeline:
    output_path: Path
    train_data: Dataset
    valid_data: Dataset
    test_data: Dataset
    trainer: Trainer

    def __init__(self, *, output_path: Path, checkpoint_path: Path | None) -> None:
        self.output_path = output_path

        NUM_STRUC = 64
        STRUC_LEN = 128
        MAX_SEQ_LEN = 32
        self.train_data = RnaDataset(
            num_struc=NUM_STRUC,
            struc_len=STRUC_LEN,
            max_seq_len=MAX_SEQ_LEN,
        )
        self.valid_data = RnaDataset(
            num_struc=NUM_STRUC,
            struc_len=STRUC_LEN,
            max_seq_len=MAX_SEQ_LEN,
        )
        self.test_data = RnaDataset(
            num_struc=NUM_STRUC,
            struc_len=STRUC_LEN,
            max_seq_len=MAX_SEQ_LEN,
        )

        model = RnaTransformer(max_seq_len=MAX_SEQ_LEN)
        cfg = TrainConfig(
            model_name="RnaTransformer",
            output_path=output_path,
            profile=False,
            save_graph=False,
            checkpoint_valid_loss=True,
            batch_size=64,
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

    def predict(self, db: str) -> None:
        primary = torch.Tensor([BOS_IDX])
        secondary = RnaTensor.from_secondary(db_to_secondary(db))
        print(secondary)
        dm = self.trainer.optimizer.dm
        # Greedy for now.
        for i in range(len(secondary)):
            out = dm([secondary.unsqueeze(0), primary.unsqueeze(0)])
            print(i)
            print(out)

    def train(self, epochs: int) -> None:
        click.echo(f"Running transformer at {self.output_path}")
        self.trainer.run(epochs)
