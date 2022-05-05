from typing import Any
from typing import Iterator
from typing import OrderedDict

from rnapy.design.harness.model import Model
import torch
from torch import nn


class DeviceModel:
    """Utilities for interacting with a model on a specified device."""

    device: str
    model: Model

    def __init__(
        self,
        *,
        model: Model,
        device: str | None = None,
    ) -> None:
        if device is None:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"
        else:
            self.device = device

        self.model = model.to(self.device)

    def __call__(self, batch: list[torch.Tensor]) -> Any:
        return self.model(*self.inputs(batch))

    def parameters(self) -> Iterator[nn.Parameter]:
        return self.model.parameters()

    def train(self) -> None:
        self.model.train()

    def eval(self) -> None:
        self.model.eval()

    def inputs(self, batch: list[torch.Tensor]) -> list[torch.Tensor]:
        return [v.to(self.device) for v in batch]

    def prediction(self, out: torch.Tensor) -> torch.Tensor:
        return self.model.model_prediction(out)

    def loss(
        self,
        *,
        batch: list[torch.Tensor],
        out: torch.Tensor,
        loss_fn: nn.Module,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        return self.model.model_loss(
            batch=self.inputs(batch), out=out.to(self.device), loss_fn=loss_fn,
        )

    def state_dict(self) -> OrderedDict[str, torch.Tensor]:
        return self.model.state_dict()

    def load_state_dict(self, state_dict: OrderedDict[str, torch.Tensor]) -> None:
        self.model.load_state_dict(state_dict)
