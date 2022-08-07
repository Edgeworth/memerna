import logging
import time

from rnapy.design.harness.config import TrainConfig
from rnapy.design.harness.trainer_protocol import TrainerProtocol
import torch
from torch import nn
from torch.utils.tensorboard import SummaryWriter


class Metrics:
    loss_sum: float = 0.0
    loss_count: float = 0.0
    accuracy_sum: float = 0.0
    accuracy_count: float = 0.0
    sample_count: float = 0.0
    time_since: float = 0.0

    def reset(self) -> None:
        self.loss_sum = 0.0
        self.loss_count = 0.0
        self.accuracy_sum = 0.0
        self.accuracy_count = 0.0
        self.sample_count = 0.0
        self.time_since = time.time()

    def record(self, loss: torch.Tensor, accuracy: torch.Tensor, sample_count: int) -> None:
        self.loss_sum += loss.sum().item()
        self.loss_count += loss.numel()
        self.accuracy_sum += accuracy.sum().item()
        self.accuracy_count += accuracy.numel()
        self.sample_count += sample_count

    def take(self) -> tuple[float, float, float]:
        loss = self.loss_sum / self.loss_count if self.loss_count != 0.0 else 0.0
        accuracy = self.accuracy_sum / self.accuracy_count if self.accuracy_count != 0.0 else 0.0
        sample_time_ms = (
            (time.time() - self.time_since) / self.sample_count if self.sample_count != 0.0 else 0.0
        )
        self.reset()
        return loss, accuracy, 1000.0 * sample_time_ms


class Reporter:
    """Reports info to tensorboard, the console, etc."""

    writer: SummaryWriter
    profiler: torch.profiler.profile | None = None
    cfg: TrainConfig

    # Tracking data:
    print_metrics: Metrics
    prev_print_sample_count: int = 0

    report_metrics: Metrics
    prev_report_sample_count: int = 0

    prev_checkpoint_sample_count: int = 0

    best_valid_loss: float = float("inf")

    # State dict elements:
    sample_count: int = 0
    sample_count_since_epoch: int = 0

    def __init__(self, cfg: TrainConfig) -> None:
        writer_path = cfg.output_path / "tensorboard" / cfg.name
        self.writer = SummaryWriter(writer_path)

        if cfg.profile:
            # Only runs one time at the beginning.
            self.profiler = torch.profiler.profile(
                schedule=torch.profiler.schedule(wait=2, warmup=2, active=6, repeat=1),
                on_trace_ready=torch.profiler.tensorboard_trace_handler(str(writer_path)),
                record_shapes=True,
                profile_memory=True,
                with_stack=False,  # TODO(1): causes segfault in pytorch 1.11.0
            )
        self.cfg = cfg

        self.print_metrics = Metrics()
        self.report_metrics = Metrics()

    def step(
        self,
        loss: torch.Tensor,
        accuracy: torch.Tensor,
        sample_count: int,
        trainer: TrainerProtocol,
    ) -> None:
        """Report metrics to tensorboard, the console, etc."""
        self.sample_count += sample_count
        self.sample_count_since_epoch += sample_count

        self.print_metrics.record(loss, accuracy, sample_count)
        self.report_metrics.record(loss, accuracy, sample_count)

        if self.sample_count - self.prev_print_sample_count >= self.cfg.print_interval:
            r_loss, r_accuracy, r_sample_time_ms = self.print_metrics.take()
            logging.info(
                f"loss: {r_loss:>7f} | accuracy: {100*r_accuracy:.2f}% | "
                f"{r_sample_time_ms:.2f} ms/sample | "
                f"{self.sample_count_since_epoch:>5d}/{self.cfg.train_samples} samples | "
                f"{self.sample_count:>5d} samples",
            )
            self.prev_print_sample_count = self.sample_count

        if self.sample_count - self.prev_report_sample_count >= self.cfg.report_interval:
            r_loss, r_accuracy, r_sample_time_ms = self.report_metrics.take()

            valid_loss, valid_accuracy = trainer.validate(self.cfg.fast_valid_samples)
            self.writer.add_scalars(
                "loss",
                {"train": r_loss, "valid": valid_loss},
                self.sample_count,
            )
            self.writer.add_scalars(
                "accuracy",
                {"train": r_accuracy * 100, "valid": valid_accuracy * 100},
                self.sample_count,
            )
            self.prev_report_sample_count = self.sample_count

        if self.sample_count - self.prev_checkpoint_sample_count >= self.cfg.checkpoint_interval:
            trainer.save_checkpoint(self.cfg.checkpoint_path(str(self.sample_count)))
            self.prev_checkpoint_sample_count = self.sample_count

        if self.profiler is not None:
            self.profiler.step()

    def on_epoch(self, loss: float, accuracy: float, trainer: TrainerProtocol) -> None:
        self.sample_count_since_epoch = 0

        if loss < self.best_valid_loss:
            self.best_valid_loss = loss
            if self.cfg.checkpoint_valid_loss:
                trainer.save_checkpoint(self.cfg.checkpoint_path("best"))

        logging.info(
            f"Epoch validation loss: {loss:>7f} accuracy: {100*accuracy:.2f}%\n",
        )

    def start(self, module: nn.Module, batch: list[torch.Tensor]) -> None:
        if self.cfg.save_graph:
            module.eval()
            self.writer.add_graph(module, batch)

        if self.profiler:
            self.profiler.__enter__()  # pylint: disable=unnecessary-dunder-call

        self.print_metrics.reset()
        self.report_metrics.reset()

    def stop(self) -> None:
        if self.profiler:
            self.profiler.__exit__(None, None, None)
        self.writer.flush()

    def state_dict(self) -> dict:
        return {"sample_count": self.sample_count}

    def load_state_dict(self, state_dict: dict) -> None:
        self.sample_count = state_dict["sample_count"]
