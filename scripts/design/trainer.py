from pathlib import Path
import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset
from torch.optim import Optimizer, Adam
from torch.utils.tensorboard import SummaryWriter
import click

from scripts.design.nn import NeuralNetwork


class Trainer:
    batch_size: int = 64
    report_interval: int = 100
    print_interval: int = 100
    device: str = "cuda" if torch.cuda.is_available() else "cpu"

    train_dataloader: DataLoader
    valid_dataloader: DataLoader
    writer: SummaryWriter
    model: nn.Module
    optimizer: Optimizer
    loss_fn: nn.CrossEntropyLoss
    global_step: int

    def __init__(self, train_data: Dataset, valid_data: Dataset, path: Path) -> None:
        # Create data loaders.
        self.train_dataloader = DataLoader(train_data, batch_size=self.batch_size)
        self.valid_dataloader = DataLoader(valid_data, batch_size=self.batch_size)
        self.writer = SummaryWriter(path / "tensorboard")

        self.model = NeuralNetwork().to(self.device)
        self.optimizer = Adam(self.model.parameters())
        self.loss_fn = nn.CrossEntropyLoss()
        self.global_step = 0

    def _next_global_step(self) -> int:
        self.global_step += 1
        return self.global_step

    def _train_epoch(self) -> None:
        self.model.train()

        avg_loss = 0.0
        for batch, (X, y) in enumerate(self.train_dataloader):
            X, y = X.to(self.device), y.to(self.device)

            # Compute prediction error
            pred = self.model(X)
            loss = self.loss_fn(pred, y)

            # Backpropagation
            self.optimizer.zero_grad()  # Zero out the gradient buffers.
            loss.backward()  # Compute gradients
            self.optimizer.step()  #  Optimize weights

            avg_loss += loss.item()

            if batch % self.print_interval == 0:
                loss, current = loss.item(), batch * len(X)
                click.echo(f"loss: {loss:>7f}  [{current:>5d}]")

            if (batch + 1) % self.report_interval == 0:
                self.writer.add_scalars(
                    "Training vs. Validation Loss",
                    {"Training": avg_loss},
                    self._next_global_step(),
                )

    def _valid_epoch(self) -> None:
        num_batches = len(self.valid_dataloader)
        self.model.eval()
        loss = 0.0
        correct = 0.0
        total = 0.0
        with torch.no_grad():  # don't calculate gradients
            for X, y in self.valid_dataloader:
                X, y = X.to(self.device), y.to(self.device)
                pred = self.model(X)
                loss += self.loss_fn(pred, y).item()
                correct += (pred.argmax(1) == y).type(torch.float).sum().item()
                total += 1.0
        loss /= num_batches
        correct /= total
        click.echo(f"Test Error: \n Correct: {correct}, Avg loss: {loss:>8f} \n")

        self.writer.add_scalars(
            "Training vs. Validation Loss",
            {"Validation": loss},
            self.global_step,
        )

    def run(self, epochs: int) -> None:
        for t in range(epochs):
            click.echo(f"Epoch {t+1}\n-------------------------------")
            self._train_epoch()
            self._valid_epoch()
        self.writer.flush()
