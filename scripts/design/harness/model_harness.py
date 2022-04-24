from typing import Any
from typing import Iterator
from typing import OrderedDict

from scripts.design.harness.model import Model
import torch
from torch import nn


class ModelHarness:
    """Utilities for interacting with a model"""

    device: str
    model: Model

    def __init__(
        self,
        *,
        model: Model,
        device: str | None = None,
    ) -> None:
        """
        Args:
            train_data: training data
            valid_data: validation data
            path: path to save model and tensorboard logs
            clip_grad_norm: clip gradient L2 norm to this value if set
        """
        if device is None:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"
        else:
            self.device = device

        self.model = model.to(self.device)

    def __call__(self, X: torch.Tensor) -> Any:
        return self.model(*self.model_input(X=X))

    def parameters(self) -> Iterator[nn.Parameter]:
        return self.model.parameters()

    def train(self) -> None:
        self.model.train()

    def eval(self) -> None:
        self.model.eval()

    def model_input(self, *, X: torch.Tensor) -> list[torch.Tensor]:
        return [t.to(self.device) for t in self.model.model_input(X=X)]

    def prediction(self, *, out: torch.Tensor) -> torch.Tensor:
        return self.model.model_prediction(out=out)

    def loss(
        self, *, out: torch.Tensor, y: torch.Tensor, loss_fn: nn.Module,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        return self.model.model_loss(out=out.to(self.device), y=y.to(self.device), loss_fn=loss_fn)

    def state_dict(self) -> OrderedDict[str, torch.Tensor]:
        return self.model.state_dict()

    def load_state_dict(self, state_dict: OrderedDict[str, torch.Tensor]) -> None:
        self.model.load_state_dict(state_dict)
