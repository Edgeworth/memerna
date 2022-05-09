import logging

from rnapy.design.harness.config import TrainConfig
from rnapy.design.harness.device_model import DeviceModel
import torch
from torch import nn


class Optimizer:
    """Encapsulates an optimizer and learning rate schedule."""

    dm: DeviceModel
    optimizer: torch.optim.Optimizer
    scheduler: torch.optim.lr_scheduler.ReduceLROnPlateau | None = None
    loss_fn: nn.CrossEntropyLoss
    cfg: TrainConfig

    def __init__(self, cfg: TrainConfig, dm: DeviceModel) -> None:
        self.dm = dm
        self.optimizer = torch.optim.SGD(self.dm.parameters(), lr=0.001)

        if cfg.optimizer == "sgd":
            self.optimizer = torch.optim.SGD(self.dm.parameters(), lr=1.0)
            self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                self.optimizer,
                factor=0.5,
                patience=5,
            )
        elif cfg.optimizer == "adam":
            self.optimizer = torch.optim.Adam(self.dm.parameters())
        elif cfg.optimizer == "amsgrad":
            self.optimizer = torch.optim.Adam(self.dm.parameters(), amsgrad=True)
        elif cfg.optimizer == "rmsprop":
            self.optimizer = torch.optim.RMSprop(self.dm.parameters())

        self.cfg = cfg

    def lr_params(self) -> list[float]:
        lrs = []
        for param_group in self.optimizer.param_groups:
            lrs.append(float(param_group["lr"]))
        return lrs

    def eval_batch(self, batch: list[torch.Tensor]) -> tuple[torch.Tensor, torch.Tensor]:
        """Evalute a single batch of data.

        Returns:
            loss: The loss tensor.
            accuracy: The accuracy tensor."""
        self.dm.eval()
        outs = self.dm(batch)
        loss, accuracy = self.dm.loss(batch=batch, outs=outs)
        return loss, accuracy

    def train_batch(self, batch: list[torch.Tensor]) -> tuple[torch.Tensor, torch.Tensor]:
        """Perform a single backprop step.

        Returns:
            loss: The loss tensor.
            accuracy: The accuracy tensor."""
        # Compute prediction error
        self.dm.train()
        outs = self.dm(batch)
        loss, accuracy = self.dm.loss(batch=batch, outs=outs)

        # Backpropagation
        self.optimizer.zero_grad()  # Zero out the gradient buffers.
        loss.backward()  # Compute gradients

        if self.cfg.clip_grad_norm:
            nn.utils.clip_grad_norm_(self.dm.parameters(), self.cfg.clip_grad_norm)

        self.optimizer.step()  #  Optimize weights
        return loss, accuracy

    def on_epoch(self, loss: float, _accuracy: float) -> None:
        if self.scheduler is not None:
            self.scheduler.step(metrics=loss)
            logging.info(f"Using learning rates: {self.lr_params()}")

    def state_dict(self) -> dict:
        d = {
            "model": self.dm.state_dict(),
            "optimizer": self.optimizer.state_dict(),
        }
        if self.scheduler is not None:
            d["scheduler"] = self.scheduler.state_dict()
        return d

    def load_state_dict(self, state_dict: dict) -> None:
        self.dm.load_state_dict(state_dict["model"])
        self.optimizer.load_state_dict(state_dict["optimizer"])
        if state_dict["scheduler"] is not None and self.scheduler is not None:
            self.scheduler.load_state_dict(state_dict["scheduler"])