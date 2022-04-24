from datetime import datetime
import itertools
import logging
import multiprocessing
from pathlib import Path

from scripts.design.harness.model import Model
from scripts.design.harness.model_harness import ModelHarness
import torch
from torch import nn
from torch.optim import Optimizer
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
from torch.utils.tensorboard import SummaryWriter


class Trainer:
    profiler: torch.profiler.profile | None = None
    name: str
    output_path: Path
    writer: SummaryWriter
    report_interval: int = 100
    """After how many minibatches to write out data to tensorboard"""

    print_interval: int = 100
    """After how many minibatches to print summary data to console"""

    checkpoint_interval: int = 10000
    """After how many minibatches to write out a checkpoint"""

    # parallelisation for data loading
    dataloader_worker_count: int = multiprocessing.cpu_count() // 4
    batch_size: int = 64
    train_dataloader: DataLoader
    valid_dataloader: DataLoader
    valid_batches: int = -1
    """How many minibatches to run validation on. -1 for all"""

    model: ModelHarness
    optimizer: Optimizer
    loss_fn: nn.CrossEntropyLoss
    clip_grad_norm: float | None = None

    step_count: int

    def __init__(
        self,
        *,
        model: Model,
        train_data: Dataset,
        valid_data: Dataset,
        output_path: Path,
        profile: bool = False,
        name: str | None = None,
        batch_size: int = 64,
        valid_batches: int = -1,
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
        self.batch_size = batch_size
        self.valid_batches = valid_batches
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

        self.model = ModelHarness(model=model)
        # TODO: Use learning schedule. Consider momentum
        # TODO: undo learning rate
        # self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.1)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.003)
        self.loss_fn = nn.CrossEntropyLoss()
        self.clip_grad_norm = clip_grad_norm
        self.step_count = 0

        if profile:
            # Only runs one time at the beginning.
            self.profiler = torch.profiler.profile(
                schedule=torch.profiler.schedule(wait=2, warmup=2, active=6, repeat=1),
                on_trace_ready=torch.profiler.tensorboard_trace_handler(str(writer_path)),
                record_shapes=True,
                profile_memory=True,
                with_stack=False,  # TODO: causes segfault in pytorch 1.11.0
            )

    def _record_graph(self) -> None:
        self.model.eval()
        X, _ = next(iter(self.train_dataloader))
        self.writer.add_graph(self.model.model, self.model.model_input(X=X))

    def _next_step(self) -> int:
        self.step_count += 1
        return self.step_count

    def _train_epoch(self) -> None:
        self.model.train()

        avg_loss = 0.0
        logging.info(f"Training epoch on {len(self.train_dataloader)} batches")
        for batch_idx, (X, y) in enumerate(self.train_dataloader, start=1):
            # Compute prediction error
            out = self.model(X)
            loss, _ = self.model.loss(out=out, y=y, loss_fn=self.loss_fn)

            # Backpropagation
            self.optimizer.zero_grad()  # Zero out the gradient buffers.
            loss.backward()  # Compute gradients

            if self.clip_grad_norm:
                nn.utils.clip_grad_norm_(self.model.parameters(), self.clip_grad_norm)

            self.optimizer.step()  #  Optimize weights
            if self.profiler:
                self.profiler.step()
            self._next_step()

            avg_loss += loss.item()

            if batch_idx % self.print_interval == 0:
                logging.info(
                    f"train loss: {loss.item():>7f} [processed {batch_idx:>5d} batches, "
                    f"{batch_idx * X.shape[0]:>5d}] inputs",
                )

            if batch_idx % self.report_interval == 0:
                avg_loss /= self.report_interval

                valid_loss, _ = self._validate()
                self.writer.add_scalars(
                    "loss",
                    {"train": avg_loss, "valid": valid_loss},
                    self.step_count,
                )
                avg_loss = 0.0

            if batch_idx % self.checkpoint_interval == 0:
                self._save_checkpoint(
                    self.output_path / self.name / f"checkpoint-{self.step_count}.pt",
                )

    def _validate(self) -> tuple[float, float]:
        total_loss = 0.0
        total_correct = 0.0
        total_examples = 0.0
        self.model.eval()

        valid_batches = min(self.valid_batches, len(self.valid_dataloader))
        if valid_batches == -1:
            valid_batches = len(self.valid_dataloader)

        logging.info(f"validating on {valid_batches} batches")
        with torch.no_grad():  # don't calculate gradients
            for X, y in itertools.islice(self.valid_dataloader, valid_batches):
                out = self.model(X)
                loss, correct = self.model.loss(out=out, y=y, loss_fn=self.loss_fn)
                total_loss += loss.item()
                total_correct += correct.sum().item()
                total_examples += correct.numel()  # Number of total predictions
        total_loss /= valid_batches
        total_correct /= total_examples
        logging.info(f"validation loss: {total_loss:>7f} accuracy: {total_correct*100:>2f}%")
        return total_loss, total_correct

    def _save_checkpoint(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        logging.info(f"Saving checkpoint to {path}")
        torch.save(
            {
                "global_step": self.step_count,
                "model_state_dict": self.model.state_dict(),
                "optimizer_state_dict": self.optimizer.state_dict(),
            },
            path,
        )

    def load_checkpoint(self, path: Path) -> None:
        checkpoint = torch.load(path)
        self.step_count = checkpoint["global_step"]
        self.model.load_state_dict(checkpoint["model_state_dict"])
        self.optimizer.load_state_dict(checkpoint["optimizer_state_dict"])

    def run(self, epochs: int) -> None:
        if self.profiler:
            self.profiler.__enter__()

        try:
            self._record_graph()
            logging.info(f"Start training for {epochs} epochs")
            for t in range(1, epochs + 1):
                logging.info(f"Epoch {t}\n-------------------------------")
                self._train_epoch()
                valid_loss, valid_correct = self._validate()
                logging.info(
                    f"Test Error: \n Correct: {valid_correct}, Avg loss: {valid_loss:>8f} \n",
                )
        finally:
            if self.profiler:
                self.profiler.__exit__(None, None, None)
            self.writer.flush()
