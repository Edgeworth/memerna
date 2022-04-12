from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple
import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset
from torch.optim import Optimizer
from torch.utils.tensorboard import SummaryWriter
import click

from scripts.design.nn import NeuralNetwork


class Trainer:
    device: str = "cuda" if torch.cuda.is_available() else "cpu"

    profiler: Optional[torch.profiler.profile] = None
    writer_path: Path
    writer: SummaryWriter
    report_interval: int = 100  # interval to write out data to tensorboard
    print_interval: int = 100  # interval to print summary data to console

    batch_size: int = 64
    train_dataloader: DataLoader
    valid_dataloader: DataLoader

    model: nn.Module
    optimizer: Optimizer
    loss_fn: nn.CrossEntropyLoss
    global_step: int

    def __init__(
        self, train_data: Dataset, valid_data: Dataset, path: Path, profile: bool = True
    ) -> None:
        # Create data loaders.
        self.train_dataloader = DataLoader(train_data, batch_size=self.batch_size)
        self.valid_dataloader = DataLoader(valid_data, batch_size=self.batch_size)

        self.writer_path = path / "tensorboard" / datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        self.writer = SummaryWriter(self.writer_path)

        self.model = NeuralNetwork().to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=1e-3)
        self.loss_fn = nn.CrossEntropyLoss()
        self.global_step = 0

        if profile:
            self.profiler = torch.profiler.profile(
                schedule=torch.profiler.schedule(wait=2, warmup=2, active=6),
                on_trace_ready=torch.profiler.tensorboard_trace_handler(str(self.writer_path)),
                record_shapes=True,
                profile_memory=True,
                with_stack=False,  # TODO: causes segfault in pytorch 1.11.0
            )

    def _record_graph(self) -> None:
        self.model.eval()
        X, _ = next(iter(self.train_dataloader))
        self.writer.add_graph(self.model, X.to(self.device))

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
            if self.profiler:
                self.profiler.step()

            avg_loss += loss.item()

            if batch % self.print_interval == 0:
                loss, current = loss.item(), batch * len(X)
                click.echo(f"loss: {loss:>7f}  [{current:>5d}]")

            if (batch + 1) % self.report_interval == 0:
                avg_loss /= self.report_interval

                valid_loss, _ = self._validate()
                self.writer.add_scalars(
                    "loss",
                    {"train": avg_loss, "valid": valid_loss},
                    self._next_global_step(),
                )
                avg_loss = 0.0

    def _validate(self) -> Tuple[float, float]:
        num_batches = len(self.valid_dataloader)
        self.model.eval()
        loss = 0.0
        acc = 0.0
        total = 0.0
        with torch.no_grad():  # don't calculate gradients
            for X, y in self.valid_dataloader:
                X, y = X.to(self.device), y.to(self.device)
                pred = self.model(X)
                loss += self.loss_fn(pred, y).item()
                acc += (pred.argmax(1) == y).type(torch.float).sum().item()
                total += 1.0
        loss /= num_batches
        acc /= total
        return loss, acc

    def run(self, epochs: int) -> None:
        if self.profiler:
            self.profiler.__enter__()

        try:
            self._record_graph()
            for t in range(epochs):
                click.echo(f"Epoch {t+1}\n-------------------------------")
                self._train_epoch()
                valid_loss, valid_acc = self._validate()
                click.echo(f"Test Error: \n Correct: {valid_acc}, Avg loss: {valid_loss:>8f} \n")
        finally:
            if self.profiler:
                self.profiler.__exit__(None, None, None)
            self.writer.flush()
