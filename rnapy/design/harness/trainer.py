import itertools
import logging
from pathlib import Path

from rnapy.design.harness.config import TrainConfig
from rnapy.design.harness.device_model import DeviceModel
from rnapy.design.harness.model import Model
from rnapy.design.harness.optimizer import Optimizer
from rnapy.design.harness.reporter import Metrics
from rnapy.design.harness.reporter import Reporter
import torch
from torch.utils.data import DataLoader
from torch.utils.data import Dataset


class Trainer:
    train_loader: DataLoader
    valid_loader: DataLoader

    optimizer: Optimizer
    reporter: Reporter
    cfg: TrainConfig

    def __init__(
        self,
        *,
        model: Model,
        train_data: Dataset,
        valid_data: Dataset,
        cfg: TrainConfig,
        checkpoint_path: Path | None = None,
    ) -> None:
        """
        Args:
            model: Model to train.
            train_data: Dataset to train on. Tensors should be of shape (batch_size, ...).
            valid_data: Dataset to validate on. Tensors should be of shape (batch_size, ...).
            cfg: Training configuration.
            checkpoint_path: Path to save checkpoint to.
        """
        # Create data loaders.
        self.train_loader = DataLoader(
            train_data,
            batch_size=cfg.batch_size,
            shuffle=True,  # Shuffle data order
            num_workers=cfg.dataloader_worker_count,  # Load data faster
        )
        self.valid_loader = DataLoader(
            valid_data,
            batch_size=cfg.batch_size,
            shuffle=True,  # Shuffle data order
            num_workers=cfg.dataloader_worker_count,  # Load data faster
        )

        checkpoint = None
        if checkpoint_path is not None:
            checkpoint = torch.load(checkpoint_path)
            cfg.load_state_dict(checkpoint["config"])

        # Set config sample counts to the actual sample counts.
        cfg.set_samples(
            len(self.train_loader) * cfg.batch_size,
            len(self.valid_loader) * cfg.batch_size,
        )

        self.optimizer = Optimizer(cfg, DeviceModel(model=model))
        self.reporter = Reporter(cfg)

        if checkpoint is not None:
            self.reporter.load_state_dict(checkpoint["reporter"])
            self.optimizer.load_state_dict(checkpoint["optimizer"])

        self.cfg = cfg

    def save_checkpoint(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        logging.info(f"Saving checkpoint to {path}")
        data = {
            "reporter": self.reporter.state_dict(),
            "optimizer": self.optimizer.state_dict(),
            "config": self.cfg.state_dict(),
        }
        torch.save(data, path)

    def _train_epoch(self) -> None:
        logging.info(f"Training epoch on {self.cfg.train_samples} samples")
        sample_count = 0
        while True:
            for batch in self.train_loader:
                loss, accuracy = self.optimizer.train_batch(batch)
                sample_count += len(batch[0])
                self.reporter.step(loss, accuracy, len(batch[0]), self)

                if sample_count > self.cfg.train_samples:
                    break
            if sample_count > self.cfg.train_samples:
                break
            logging.warning(f"reached end of train data after {sample_count} samples. reusing...")

    def validate(self, num_samples: int) -> tuple[float, float]:
        logging.info(f"validating on {num_samples} samples")
        metrics = Metrics()
        sample_count = 0
        with torch.no_grad():  # don't calculate gradients
            while True:
                for batch in self.valid_loader:
                    loss, accuracy = self.optimizer.eval_batch(batch)
                    sample_count += len(batch[0])
                    metrics.record(loss, accuracy, len(batch[0]))

                    if sample_count > num_samples:
                        break
                if sample_count > num_samples:
                    break
                logging.warning(
                    f"reached end of validation data after {sample_count} samples. reusing..."
                )
        r_loss, r_accuracy, _ = metrics.take()
        logging.info(f"validation loss: {r_loss:>7f} accuracy: {r_accuracy*100:.2f}%")
        return r_loss, r_accuracy

    def run(self, epochs: int) -> None:
        # Add graph of model to tensorboard.
        dm = self.optimizer.dm
        dm.eval()
        batch = next(iter(self.train_loader))
        self.reporter.start(dm.model, dm.inputs(batch))

        try:
            logging.info(
                f"Start training for {epochs} epochs, model has {dm.num_parameters()} parameters"
            )
            for t in range(1, epochs + 1):
                logging.info(f"Epoch {t}\n-------------------------------")
                self._train_epoch()
                loss, accuracy = self.validate(self.cfg.accurate_valid_samples)
                self.reporter.on_epoch(loss, accuracy, self)
                self.optimizer.on_epoch(loss, accuracy)

        finally:
            self.reporter.stop()
