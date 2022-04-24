import logging

from scripts.design.harness.config import TrainConfig
import torch
from torch import nn
from torch.utils.tensorboard import SummaryWriter

from scripts.design.harness.trainer_protocol import TrainerProtocol


class Reporter:
    """Reports info to tensorboard, the console, etc."""

    writer: SummaryWriter
    profiler: torch.profiler.profile | None = None
    cfg: TrainConfig

    # State dict elements:
    step_count: int = 1
    loss_sum: float = 0.0
    loss_count: float = 0.0
    accuracy_sum: float = 0.0
    accuracy_count: float = 0.0

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

    def step(
        self,
        loss: torch.Tensor,
        accuracy: torch.Tensor,
        trainer: TrainerProtocol,
        batch_idx: int,
    ) -> None:
        """Report metrics to tensorboard, the console, etc."""
        self.loss_sum += loss.sum().item()
        self.loss_count += loss.numel()
        self.accuracy_sum += accuracy.sum().item()
        self.accuracy_count += accuracy.numel()

        if self.step_count % self.cfg.print_interval == 0:
            logging.info(
                f"train loss: {loss.item():>7f} [processed "
                f"{batch_idx:>5d}/{self.cfg.train_batches} batches]",
            )

        if self.step_count % self.cfg.report_interval == 0:
            self.loss_sum /= self.loss_count
            self.accuracy_sum /= self.accuracy_count

            valid_loss, valid_accuracy = trainer.validate(self.cfg.fast_valid_batches)
            self.writer.add_scalars(
                "loss",
                {"train": self.loss_sum, "valid": valid_loss},
                self.step_count,
            )
            self.writer.add_scalars(
                "accuracy",
                {"train": self.accuracy_sum * 100, "valid": valid_accuracy * 100},
                self.step_count,
            )
            self.loss_sum = 0.0
            self.loss_count = 0.0
            self.accuracy_sum = 0.0
            self.accuracy_count = 0.0

        if self.step_count % self.cfg.checkpoint_interval == 0:
            trainer.save_checkpoint(self.cfg.checkpoint_path(self.step_count))

        if self.profiler is not None:
            self.profiler.step()

        self.step_count += 1

    def on_epoch(self, loss: float, accuracy: float) -> None:
        logging.info(
            f"Epoch validation loss: {loss:>7f} accuracy: {100*accuracy:.2f}%\n",
        )

    def start(self, module: nn.Module, inputs: list[torch.Tensor]) -> None:
        module.eval()
        self.writer.add_graph(module, inputs)

        if self.profiler:
            self.profiler.__enter__()

    def stop(self) -> None:
        if self.profiler:
            self.profiler.__exit__(None, None, None)
        self.writer.flush()

    def state_dict(self) -> dict:
        return {
            "step_count": self.step_count,
            "avg_loss": self.loss_sum,
        }

    def load_state_dict(self, state_dict: dict) -> None:
        self.step_count = state_dict["step_count"]
        self.loss_sum = state_dict["avg_loss"]
