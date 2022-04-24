from datetime import datetime
import itertools
import logging
import multiprocessing
from pathlib import Path

from scripts.design.harness.model import Model
from scripts.design.harness.model_harness import ModelHarness
import torch
from torch import nn
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
    train_batches: int = -1
    """How many minibatches to run training on per epoch. -1 for all"""

    valid_dataloader: DataLoader
    fast_valid_batches: int = -1
    """How many minibatches to run validation on for a quick guess (e.g.
    reporting validation loss regularly). -1 for all"""

    accurate_valid_batches: int = -1
    """How many minibatches to run validation on for accuracy (e.g. at end of
    each epoch). -1 for all"""

    model: ModelHarness
    optimizer: torch.optim.Optimizer
    scheduler: torch.optim.lr_scheduler.ReduceLROnPlateau | None
    loss_fn: nn.CrossEntropyLoss
    clip_grad_norm: float | None = None

    step_count: int = 0

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
        train_batches: int = -1,
        fast_valid_batches: int = -1,
        accurate_valid_batches: int = -1,
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
        self.output_path = output_path

        writer_path = output_path / "tensorboard" / self.name
        self.writer = SummaryWriter(writer_path)

        # Create data loaders.
        def batch_num(batches: int, loader: DataLoader) -> int:
            if batches == -1:
                return len(loader)
            return min(batches, len(loader))

        self.train_dataloader = DataLoader(
            train_data,
            batch_size=self.batch_size,
            shuffle=True,  # Shuffle data order
            num_workers=self.dataloader_worker_count,  # Load data faster
        )
        self.train_batches = batch_num(train_batches, self.train_dataloader)

        self.valid_dataloader = DataLoader(
            valid_data,
            batch_size=self.batch_size,
            shuffle=True,  # Shuffle data order
            num_workers=self.dataloader_worker_count,  # Load data faster
        )
        self.fast_valid_batches = batch_num(fast_valid_batches, self.valid_dataloader)
        self.accurate_valid_batches = batch_num(accurate_valid_batches, self.valid_dataloader)

        self.model = ModelHarness(model=model)
        # TODO: Use learning schedule. Consider momentum
        # TODO: undo learning rate
        self.optimizer = torch.optim.SGD(self.model.parameters(), lr=1.0)
        # self.optimizer = torch.optim.Adam(self.model.parameters())

        if True:
            self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                self.optimizer, factor=0.1, patience=2
            )

        self.loss_fn = nn.CrossEntropyLoss()
        self.clip_grad_norm = clip_grad_norm

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

    def _train_batch(self, X: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
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

        return loss

    def _train_epoch(self) -> None:
        self.model.train()

        avg_loss = 0.0
        logging.info(f"Training epoch on {self.train_batches} batches")
        for batch_idx, (X, y) in enumerate(
            itertools.islice(self.train_dataloader, self.train_batches), start=1
        ):
            loss = self._train_batch(X, y)

            avg_loss += loss.item()

            if self.step_count % self.print_interval == 0:
                logging.info(
                    f"train loss: {loss.item():>7f} [processed "
                    f"{batch_idx:>5d}/{self.train_batches} batches, "
                    f"{batch_idx * X.shape[0]:>7d} inputs, {self.step_count:>7d} total batches]",
                )

            if self.step_count % self.report_interval == 0:
                avg_loss /= self.report_interval

                valid_loss, _ = self._validate(self.fast_valid_batches)
                self.writer.add_scalars(
                    "loss",
                    {"train": avg_loss, "valid": valid_loss},
                    self.step_count,
                )
                avg_loss = 0.0

            if self.step_count % self.checkpoint_interval == 0:
                self._save_checkpoint(
                    self.output_path / self.name / f"checkpoint-{self.step_count}.pt",
                )

    def _validate(self, num_batches: int) -> tuple[float, float]:
        total_loss = 0.0
        total_correct = 0.0
        total_examples = 0.0
        self.model.eval()
        logging.info(f"validating on {num_batches} batches")
        with torch.no_grad():  # don't calculate gradients
            for X, y in itertools.islice(self.valid_dataloader, num_batches):
                out = self.model(X)
                loss, correct = self.model.loss(out=out, y=y, loss_fn=self.loss_fn)
                total_loss += loss.item()
                total_correct += correct.sum().item()
                total_examples += correct.numel()  # Number of total predictions
        total_loss /= num_batches
        total_correct /= total_examples
        logging.info(f"validation loss: {total_loss:>7f} accuracy: {total_correct*100:.2f}%")
        return total_loss, total_correct

    def _save_checkpoint(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        logging.info(f"Saving checkpoint to {path}")
        data = {
            "global_step": self.step_count,
            "model_state_dict": self.model.state_dict(),
            "optimizer_state_dict": self.optimizer.state_dict(),
        }
        if self.scheduler is not None:
            data["scheduler_state_dict"] = self.scheduler.state_dict()
        torch.save(data, path)

    def load_checkpoint(self, path: Path) -> None:
        checkpoint = torch.load(path)
        self.step_count = checkpoint["global_step"]
        self.model.load_state_dict(checkpoint["model_state_dict"])
        self.optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
        if "scheduler_state_dict" in checkpoint and self.scheduler is not None:
            self.scheduler.load_state_dict(checkpoint["scheduler_state_dict"])

    def _learning_rate_parameters(self) -> list[float]:
        lrs = []
        for param_group in self.optimizer.param_groups:
            lrs.append(float(param_group["lr"]))
        return lrs

    def run(self, epochs: int) -> None:
        if self.profiler:
            self.profiler.__enter__()

        try:
            self._record_graph()
            logging.info(f"Start training for {epochs} epochs")
            for t in range(1, epochs + 1):
                logging.info(f"Epoch {t}\n-------------------------------")
                self._train_epoch()
                valid_loss, valid_correct = self._validate(self.accurate_valid_batches)
                if self.scheduler is not None:
                    self.scheduler.step(metrics=valid_loss)
                    logging.info(f"Using learning rates: {self._learning_rate_parameters()}")

                logging.info(
                    f"Epoch validation loss: {valid_loss:>7f} accuracy: {100*valid_correct:.2f}%\n",
                )
        finally:
            if self.profiler:
                self.profiler.__exit__(None, None, None)
            self.writer.flush()
