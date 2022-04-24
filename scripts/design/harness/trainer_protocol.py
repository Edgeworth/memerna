from pathlib import Path
from typing import Protocol


class TrainerProtocol(Protocol):
    def save_checkpoint(self, path: Path) -> None:
        raise NotImplementedError

    def validate(self, num_batches: int) -> tuple[float, float]:
        raise NotImplementedError
