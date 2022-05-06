from dataclasses import dataclass
from dataclasses import field
from datetime import datetime
import multiprocessing
from pathlib import Path


@dataclass(eq=True, kw_only=True)
class TrainConfig:
    """Encapsulates all configuration for running training for a model."""

    model_name: str
    """Name of the model"""

    output_path: Path
    """Path to save model and tensorboard logs"""

    batch_size: int = 64
    """Batch size for training"""

    train_samples: int = -1
    """How many samples to run training on per epoch. -1 for all"""

    fast_valid_samples: int = -1
    """How many samples to run validation on for a quick guess (e.g.
    reporting validation loss regularly). -1 for all"""

    accurate_valid_samples: int = -1
    """How many samples to run validation on for accuracy (e.g. at end of
    each epoch). -1 for all"""

    dataloader_worker_count: int = multiprocessing.cpu_count() // 4
    """How many workers to use for data loading."""

    report_interval: int = 1024
    """After how many samples to write out data to tensorboard"""

    print_interval: int = 1024
    """After how many samples to print summary data to console"""

    checkpoint_interval: int = 200000
    """After how many samples to write out a checkpoint"""

    checkpoint_valid_loss: bool = False
    """Whether to checkpoint the best validation loss model so far."""

    save_graph: bool = False
    """Whether to save the model's graph to tensorboard."""

    profile: bool = False
    """Whether to enable profiling of the model"""

    optimizer: str = "sgd"
    """Which optimizer to use"""

    # TODO: Try autocast.

    clip_grad_norm: float | None = None
    """Clip gradient L2 norm to this value if set"""

    name: str = field(init=False)
    """Automatically generated name of the config"""

    def __post_init__(self) -> None:
        # Don't include anything except the model name and time, to allow
        # a model to be trained with multiple configs (using checkpointing).
        time_str = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        self.name = f"{self.model_name}-{time_str}"

    def set_samples(self, max_train_samples: int, max_valid_samples: int) -> None:
        def sample_num(samples: int, max_samples: int) -> int:
            if samples == -1:
                return max_samples
            return min(samples, max_samples)

        self.train_samples = sample_num(self.train_samples, max_train_samples)
        self.fast_valid_samples = sample_num(self.fast_valid_samples, max_valid_samples)
        self.accurate_valid_samples = sample_num(self.accurate_valid_samples, max_valid_samples)

    def checkpoint_path(self, tag: str) -> Path:
        return self.output_path / self.name / f"checkpoint-{tag}.pt"

    def state_dict(self) -> dict:
        # Just save the config name, so we can restart from the same e.g. tensorboard log dir.
        return {"name": self.name}

    def load_state_dict(self, state_dict: dict) -> None:
        self.name = state_dict["name"]
