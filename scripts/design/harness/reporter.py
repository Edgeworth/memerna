import logging
import time

from scripts.design.harness.config import TrainConfig
import torch
from torch import nn
from torch.utils.tensorboard import SummaryWriter

from scripts.design.harness.trainer_protocol import TrainerProtocol


class Metrics:
    loss_sum: float = 0.0
    loss_count: float = 0.0
    accuracy_sum: float = 0.0
    accuracy_count: float = 0.0
    record_count: float = 0.0
    time_since: float = 0.0

    def reset(self) -> None:
        self.loss_sum = 0.0
        self.loss_count = 0.0
        self.accuracy_sum = 0.0
        self.accuracy_count = 0.0
        self.record_count = 0.0
        self.time_since = time.time()

    def record(self, loss: torch.Tensor, accuracy: torch.Tensor) -> None:
        self.loss_sum += loss.sum().item()
        self.loss_count += loss.numel()
        self.accuracy_sum += accuracy.sum().item()
        self.accuracy_count += accuracy.numel()
        self.record_count += 1

    def take(self, cfg: TrainConfig) -> tuple[float, float, float]:
        loss = self.loss_sum / self.loss_count if self.loss_count != 0.0 else 0.0
        accuracy = self.accuracy_sum / self.accuracy_count if self.accuracy_count != 0.0 else 0.0
        batch_time_ms = (
            (time.time() - self.time_since) / self.record_count if self.record_count != 0.0 else 0.0
        )
        batch_time_ms = 1000.0 * batch_time_ms / cfg.batch_size
        self.reset()
        return loss, accuracy, batch_time_ms


class Reporter:
    """Reports info to tensorboard, the console, etc."""

    writer: SummaryWriter
    profiler: torch.profiler.profile | None = None
    cfg: TrainConfig

    # Tracking data:
    print_metrics: Metrics
    report_metrics: Metrics

    # State dict elements:
    step_count: int = 0
    step_count_since_epoch: int = 0

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
                with_stack=False,  # TODO: causes segfault in pytorch 1.11.0
            )
        self.cfg = cfg

        self.print_metrics = Metrics()
        self.report_metrics = Metrics()

    def step(
        self,
        loss: torch.Tensor,
        accuracy: torch.Tensor,
        trainer: TrainerProtocol,
    ) -> None:
        """Report metrics to tensorboard, the console, etc."""
        self.step_count += 1
        self.step_count_since_epoch += 1

        self.print_metrics.record(loss, accuracy)
        self.report_metrics.record(loss, accuracy)

        if self.step_count % self.cfg.print_interval == 0:
            r_loss, r_accuracy, r_batch_time_ms = self.print_metrics.take(self.cfg)
            logging.info(
                f"train loss: {r_loss:>7f} | train accuracy: {100*r_accuracy:.2f}% | "
                f"{r_batch_time_ms:.2f} ms/batch | "
                f"{self.step_count_since_epoch:>5d}/{self.cfg.train_batches} batches | "
                f"{self.step_count:>5d} steps",
            )

        if self.step_count % self.cfg.report_interval == 0:
            r_loss, r_accuracy, r_batch_time_ms = self.report_metrics.take(self.cfg)

            valid_loss, valid_accuracy = trainer.validate(self.cfg.fast_valid_batches)
            self.writer.add_scalars(
                "loss",
                {"train": r_loss, "valid": valid_loss},
                self.step_count,
            )
            self.writer.add_scalars(
                "accuracy",
                {"train": r_accuracy * 100, "valid": valid_accuracy * 100},
                self.step_count,
            )
            self.writer.add_scalars("batch time ms", {"train": r_batch_time_ms}, self.step_count)

        if self.step_count % self.cfg.checkpoint_interval == 0:
            trainer.save_checkpoint(self.cfg.checkpoint_path(self.step_count))

        if self.profiler is not None:
            self.profiler.step()

    def on_epoch(self, loss: float, accuracy: float) -> None:
        self.step_count_since_epoch = 0
        logging.info(
            f"Epoch validation loss: {loss:>7f} accuracy: {100*accuracy:.2f}%\n",
        )

    def start(self, module: nn.Module, inputs: list[torch.Tensor]) -> None:
        if self.cfg.save_graph:
            module.eval()
            self.writer.add_graph(module, inputs)

        if self.profiler:
            self.profiler.__enter__()

        self.print_metrics.reset()
        self.report_metrics.reset()

    def stop(self) -> None:
        if self.profiler:
            self.profiler.__exit__(None, None, None)
        self.writer.flush()

    def state_dict(self) -> dict:
        return {
            "step_count": self.step_count,
        }

    def load_state_dict(self, state_dict: dict) -> None:
        self.step_count = state_dict["step_count"]
