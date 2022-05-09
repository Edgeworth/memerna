from pathlib import Path

import click
from rnapy.design.harness.config import TrainConfig
from rnapy.design.harness.trainer import Trainer
from rnapy.design.rna.config import RnaPipelineConfig
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

    def __init__(
        self, *, cfg: RnaPipelineConfig, output_path: Path, checkpoint_path: Path | None
    ) -> None:
        self.output_path = output_path

        self.train_data = RnaDataset(num_struc=cfg.train_num_struc, cfg=cfg)
        self.valid_data = RnaDataset(num_struc=cfg.valid_num_struc, cfg=cfg)
        self.test_data = RnaDataset(num_struc=cfg.valid_num_struc, cfg=cfg)

        model = RnaTransformer(cfg=cfg)
        train_cfg = TrainConfig(
            model_name="RnaTransformer",
            output_path=output_path,
            profile=False,
            save_graph=False,
            checkpoint_valid_loss=True,
            optimizer="adam",
            batch_size=cfg.batch_size,
            train_samples=10000,
            fast_valid_samples=512,
            accurate_valid_samples=4096,
            clip_grad_norm=1.0,
        )
        self.trainer = Trainer(
            model=model,
            train_data=self.train_data,
            valid_data=self.valid_data,
            cfg=train_cfg,
            checkpoint_path=checkpoint_path,
        )

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
