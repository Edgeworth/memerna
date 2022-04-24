from typing import Any

from scripts.design.harness.model import Model
import torch
from torch import nn


class FashionModel(Model):
    def __init__(self) -> None:
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(28 * 28, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, 10),
        )

    def forward(self, x: torch.Tensor) -> Any:
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits
