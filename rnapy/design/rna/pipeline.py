from pathlib import Path

import click
from rnapy.design.harness.train_cfg import TrainCfg
from rnapy.design.harness.trainer import Trainer
from rnapy.design.rna.dataset import RnaDataset
from rnapy.design.rna.pipeline_cfg import RnaPipelineCfg
from rnapy.design.rna.tensor import BOS_IDX
import torch
from torch.utils.data import Dataset


class RnaPipeline:
    output_path: Path
    train_data: Dataset
    valid_data: Dataset
    test_data: Dataset
    trainer: Trainer
    cfg: RnaPipelineCfg

    def __init__(
        self,
        *,
        cfg: RnaPipelineCfg,
        output_path: Path,
        checkpoint_path: Path | None,
    ) -> None:
        self.output_path = output_path
        self.cfg = cfg

        # TODO(2): Inference without all this training infrastructure.
        self.train_data = RnaDataset(num_struc=cfg.train_num_struc, cfg=cfg)
        self.valid_data = RnaDataset(num_struc=cfg.valid_num_struc, cfg=cfg)
        self.test_data = RnaDataset(num_struc=cfg.valid_num_struc, cfg=cfg)

        model = cfg.model_class(cfg=cfg)
        # TODO(1): Defaults for train config per model kind? e.g. checkpoint_interval
        train_cfg = TrainCfg(
            model_name="RnaPipeline",
            output_path=output_path,
            profile=False,
            save_graph=False,
            checkpoint_valid_loss=True,
            optimizer="adam",
            batch_size=cfg.batch_size,
            train_samples=cfg.batch_size * 128,  # 128 batches per epoch.
            fast_valid_samples=cfg.batch_size * 4,  # 4 batches
            accurate_valid_samples=cfg.batch_size * 32,  # 32 batches
            clip_grad_norm=1.0,
        )
        self.trainer = Trainer(
            model=model,
            train_data=self.train_data,
            valid_data=self.valid_data,
            cfg=train_cfg,
            checkpoint_path=checkpoint_path,
        )

    def predict_db_only(self, db: str) -> str:
        """Predicts a primary from a dot-bracket sequence only."""
        db_tensor = self.cfg.tensor.from_db(db)
        dm = self.trainer.optimizer.dm
        out = dm([db_tensor.unsqueeze(0), torch.empty()])
        pred = dm.prediction(out=out[0]).to("cpu").squeeze(0)
        return self.cfg.tensor.to_primary(pred)

    def predict_mlm(self, db: str) -> str:
        primary = [BOS_IDX]
        db_tensor = self.cfg.tensor.from_db(db)
        print(db_tensor)
        dm = self.trainer.optimizer.dm
        # Greedy for now.
        for _ in range(len(db_tensor)):
            out = dm([db_tensor.unsqueeze(0), torch.Tensor(primary).unsqueeze(0)])
            pred = dm.prediction(out=out).to("cpu").squeeze(0)
            print(primary, pred[-1].item())
            primary.append(int(pred[-1].item()))
        return self.cfg.tensor.to_primary(primary)

    def train(self, epochs: int) -> None:
        click.echo(f"Running transformer at {self.output_path}")
        self.trainer.run(epochs)
