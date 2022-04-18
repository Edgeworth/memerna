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

from scripts.design.harness.model import Model


class Trainer:
    device: str = "cuda" if torch.cuda.is_available() else "cpu"

    profiler: torch.profiler.profile | None = None
    name: str
    output_path: Path
    writer: SummaryWriter
    report_interval: int = 100  # interval to write out data to tensorboard
    print_interval: int = 100  # interval to print summary data to console

    # parallelisation for data loading
    dataloader_worker_count: int = multiprocessing.cpu_count() // 4
    batch_size: int = 64
    train_dataloader: DataLoader
    valid_dataloader: DataLoader

    model: Model
    optimizer: Optimizer
    loss_fn: nn.CrossEntropyLoss
    clip_grad_norm: float | None = None
    global_step: int

    def __init__(
        self,
        *,
        model: Model,
        train_data: Dataset,
        valid_data: Dataset,
        output_path: Path,
        name: str | None = None,
        profile: bool = False,
        clip_grad_norm: float | None = None,
    ) -> None:
        """
        Args:
            train_data: training data
            valid_data: validation data
            path: path to save model and tensorboard logs
            clip_grad_norm: clip gradient L2 norm to this value if set
        """
        time_str = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        self.name = name or f"{model.__class__.__name__}-{time_str}"
        self.output_path = output_path

        writer_path = output_path / "tensorboard" / self.name
        self.writer = SummaryWriter(writer_path)

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
                on_trace_ready=torch.profiler.tensorboard_trace_handler(str(writer_path)),
                record_shapes=True,
                profile_memory=False,
                with_stack=False,  # TODO: causes segfault in pytorch 1.11.0
            )

    def _model_input(self, *, X: torch.Tensor) -> list[torch.Tensor]:
        return [t.to(self.device) for t in self.model.model_input(X=X)]

    def _record_graph(self) -> None:
        self.model.eval()
        X, _ = next(iter(self.train_dataloader))
        self.writer.add_graph(self.model, self._model_input(X=X))

    def _next_global_step(self) -> int:
        self.global_step += 1
        return self.global_step

    def _train_epoch(self) -> None:
        self.model.train()

        avg_loss = 0.0
        for batch_idx, (X, y) in enumerate(self.train_dataloader, start=1):
            X, y = X.to(self.device), y.to(self.device)

            # Compute prediction error
            out = self.model(*self._model_input(X=X))
            loss, _ = self.model.model_loss(out=out, y=y, loss_fn=self.loss_fn)

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
                logging.info(f"loss: {loss.item():>7f}  [processed {batch_idx * len(X):>5d}]")

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
        total_loss = 0.0
        total_correct = 0.0
        total_examples = 0.0
        self.model.eval()
        with torch.no_grad():  # don't calculate gradients
            for X, y in self.valid_dataloader:
                X, y = X.to(self.device), y.to(self.device)
                out = self.model(*self._model_input(X=X))
                loss, correct = self.model.model_loss(out=out, y=y, loss_fn=self.loss_fn)
                total_loss += loss.item()
                total_correct += correct.sum().item()
                total_examples += len(X)
        total_loss /= len(self.valid_dataloader)
        total_correct /= total_examples
        return total_loss, total_correct

    def _save_checkpoint(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        torch.save(
            {
                "global_step": self.global_step,
                "model_state_dict": self.model.state_dict(),
                "optimizer_state_dict": self.optimizer.state_dict(),
            },
            path,
        )

    def _load_checkpoint(self, path: Path) -> None:
        checkpoint = torch.load(path)
        self.global_step = checkpoint["global_step"]
        self.model.load_state_dict(checkpoint["model_state_dict"])
        self.optimizer.load_state_dict(checkpoint["optimizer_state_dict"])

    def run(self, epochs: int) -> None:
        if self.profiler:
            self.profiler.__enter__()

        try:
            self._record_graph()
            for t in range(1, epochs + 1):
                logging.info(f"Epoch {t}\n-------------------------------")
                self._train_epoch()
                valid_loss, valid_correct = self._validate()
                self._save_checkpoint(self.output_path / self.name / f"checkpoint-{t}.pt")
                logging.info(
                    f"Test Error: \n Correct: {valid_correct}, Avg loss: {valid_loss:>8f} \n"
                )
        finally:
            if self.profiler:
                self.profiler.__exit__(None, None, None)
            self.writer.flush()
