import itertools
import logging
from pathlib import Path

from scripts.design.harness.config import TrainConfig
from scripts.design.harness.device_model import DeviceModel
from scripts.design.harness.model import Model
from scripts.design.harness.optimizer import Optimizer
from scripts.design.harness.reporter import Reporter
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

        # Set config batch counts to the actual batch counts.
        cfg.set_max_batches(len(self.train_loader), len(self.valid_loader))

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
        logging.info(f"Training epoch on {self.cfg.train_batches} batches")
        for X, y in itertools.islice(self.train_loader, self.cfg.train_batches):
            loss, accuracy = self.optimizer.train_batch(X, y)
            self.reporter.step(loss, accuracy, self)

    def validate(self, num_batches: int) -> tuple[float, float]:
        avg_loss = 0.0
        avg_accuracy = 0.0
        num_inputs = 0.0
        logging.info(f"validating on {num_batches} batches")
        with torch.no_grad():  # don't calculate gradients
            for X, y in itertools.islice(self.valid_loader, num_batches):
                loss, correct = self.optimizer.eval_batch(X, y)
                avg_loss += loss.item()
                avg_accuracy += correct.sum().item()
                num_inputs += correct.numel()  # Number of total predictions
        avg_loss /= num_batches
        avg_accuracy /= num_inputs
        logging.info(f"validation loss: {avg_loss:>7f} accuracy: {avg_accuracy*100:.2f}%")
        return avg_loss, avg_accuracy

    def run(self, epochs: int) -> None:
        # Add graph of model to tensorboard.
        dm = self.optimizer.dm
        dm.eval()
        X, _ = next(iter(self.train_loader))
        self.reporter.start(dm.model, dm.inputs(X))

        try:
            logging.info(f"Start training for {epochs} epochs")
            for t in range(1, epochs + 1):
                logging.info(f"Epoch {t}\n-------------------------------")
                self._train_epoch()
                loss, accuracy = self.validate(self.cfg.accurate_valid_batches)
                self.reporter.on_epoch(loss, accuracy)
                self.optimizer.on_epoch(loss, accuracy)

        finally:
            self.reporter.stop()
