from pathlib import Path

import click
from rnapy.design.harness.config import TrainConfig
from rnapy.design.harness.trainer import Trainer
from rnapy.design.rna.dataset import RnaDataset
from rnapy.design.rna.rna_transformer import RnaTransformer
from torch.utils.data import Dataset
from rnapy.design.rna.tensor import BOS_IDX, RnaTensor
import torch


class RnaPipeline:
    output_path: Path
    train_data: Dataset
    valid_data: Dataset
    test_data: Dataset
    trainer: Trainer

    def __init__(self, *, output_path: Path, checkpoint_path: Path | None) -> None:
        self.output_path = output_path

        STRUC_LEN = 128
        MAX_SEQ_LEN = 32
        BATCH_SIZE = 128
        self.train_data = RnaDataset(
            num_struc=1024,
            struc_len=STRUC_LEN,
            max_seq_len=MAX_SEQ_LEN,
        )
        self.valid_data = RnaDataset(
            num_struc=64,
            struc_len=STRUC_LEN,
            max_seq_len=MAX_SEQ_LEN,
        )
        self.test_data = RnaDataset(
            num_struc=64,
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
            batch_size=BATCH_SIZE,
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
        primary = [BOS_IDX]
        db_tensor = RnaTensor.from_db(db)
        print(db_tensor)
        dm = self.trainer.optimizer.dm
        # Greedy for now.
        for _ in range(len(db_tensor)):
            out = dm([db_tensor.unsqueeze(0), torch.Tensor(primary).unsqueeze(0)])
            pred = dm.prediction(out=out).to("cpu").squeeze(0)
            print(primary, pred[-1].item())
            primary.append(int(pred[-1].item()))
        print(RnaTensor.to_primary(primary))

    def train(self, epochs: int) -> None:
        click.echo(f"Running transformer at {self.output_path}")
        self.trainer.run(epochs)
