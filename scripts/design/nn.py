from typing import Any
from torch import nn

# Define model
class NeuralNetwork(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(28 * 28, 512), nn.ReLU(), nn.Linear(512, 512), nn.ReLU(), nn.Linear(512, 10)
        )

    def forward(self, x: Any) -> Any:
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits
