from datetime import datetime
import logging
import multiprocessing
from pathlib import Path

import torch
from torch import nn
from torch.optim import Optimizer
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
from torch.utils.tensorboard import SummaryWriter


class Trainer:
    device: str = "cuda" if torch.cuda.is_available() else "cpu"

    profiler: torch.profiler.profile | None = None
    writer_path: Path
    writer: SummaryWriter
    report_interval: int = 1  # interval to write out data to tensorboard
    print_interval: int = 1  # interval to print summary data to console

    # parallelisation for data loading
    dataloader_worker_count: int = multiprocessing.cpu_count() // 4
    batch_size: int = 64
    train_dataloader: DataLoader
    valid_dataloader: DataLoader

    model: nn.Module
    optimizer: Optimizer
    loss_fn: nn.CrossEntropyLoss
    clip_grad_norm: float | None = None
    global_step: int

    def __init__(
        self,
        *,
        model: nn.Module,
        train_data: Dataset,
        valid_data: Dataset,
        path: Path,
        clip_grad_norm: float | None = None,
        profile: bool = False,
    ) -> None:
        """
        Args:
            train_data: training data
            valid_data: validation data
            path: path to save model and tensorboard logs
            clip_grad_norm: clip gradient L2 norm to this value if set
        """
        # Create data loaders.
        self.train_dataloader = DataLoader(
            train_data,
            batch_size=self.batch_size,
            shuffle=True,  # Shuffle data order
            num_workers=self.dataloader_worker_count,  # Load data faster
        )
        self.valid_dataloader = DataLoader(
            valid_data,
            batch_size=self.batch_size,
            shuffle=True,  # Shuffle data order
            num_workers=self.dataloader_worker_count,  # Load data faster
        )

        self.writer_path = path / "tensorboard" / datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        self.writer = SummaryWriter(self.writer_path)

        self.model = model.to(self.device)
        # TODO: Use learning schedule. Consider momentum
        self.optimizer = torch.optim.SGD(self.model.parameters(), lr=1e-3)
        self.loss_fn = nn.CrossEntropyLoss()
        self.clip_grad_norm = clip_grad_norm
        self.global_step = 0

        if profile:
            # Only runs one time at the beginning.
            self.profiler = torch.profiler.profile(
                schedule=torch.profiler.schedule(wait=2, warmup=2, active=6, repeat=1),
                on_trace_ready=torch.profiler.tensorboard_trace_handler(str(self.writer_path)),
                record_shapes=True,
                profile_memory=False,
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
        for batch_idx, (X, y) in enumerate(self.train_dataloader, start=1):
            X, y = X.to(self.device), y.to(self.device)

            # Compute prediction error
            pred = self.model(X, X).reshape(-1, 28782)  # TODO: undo
            loss = self.loss_fn(pred, y.reshape(-1))  # TODO: undo

            # Backpropagation
            self.optimizer.zero_grad()  # Zero out the gradient buffers.
            loss.backward()  # Compute gradients

            if self.clip_grad_norm:
                nn.utils.clip_grad_norm_(self.model.parameters(), self.clip_grad_norm)

            self.optimizer.step()  #  Optimize weights
            if self.profiler:
                self.profiler.step()

            avg_loss += loss.item()

            if batch_idx % self.print_interval == 0:
                loss, current = loss.item(), batch_idx * len(X)
                logging.info(f"loss: {loss:>7f}  [{current:>5d}]")

            if batch_idx % self.report_interval == 0:
                avg_loss /= self.report_interval

                valid_loss, _ = self._validate()
                self.writer.add_scalars(
                    "loss",
                    {"train": avg_loss, "valid": valid_loss},
                    self._next_global_step(),
                )
                avg_loss = 0.0

    def _validate(self) -> tuple[float, float]:
        num_batches = len(self.valid_dataloader)
        loss = 0.0
        acc = 0.0
        total_examples = 0.0  #
        self.model.eval()
        with torch.no_grad():  # don't calculate gradients
            for X, y in self.valid_dataloader:
                X, y = X.to(self.device), y.to(self.device)
                pred = self.model(X, X).reshape(-1, 28782)  # TODO: undo
                loss += self.loss_fn(pred, y.reshape(-1)).item()  # TODO: undo
                # TODO: undo
                acc += (pred.argmax(1) == y.reshape(-1)).type(torch.float).sum().item()
                total_examples += len(X)
                break  # TODO: undo
        # TODO: undo
        # loss /= num_batches
        acc /= total_examples
        return loss, acc

    def run(self, epochs: int) -> None:
        if self.profiler:
            self.profiler.__enter__()

        try:
            # TODO: undo
            # self._record_graph()
            for t in range(epochs):
                logging.info(f"Epoch {t+1}\n-------------------------------")
                self._train_epoch()
                valid_loss, valid_acc = self._validate()
                logging.info(f"Test Error: \n Correct: {valid_acc}, Avg loss: {valid_loss:>8f} \n")
        finally:
            if self.profiler:
                self.profiler.__exit__(None, None, None)
            self.writer.flush()
