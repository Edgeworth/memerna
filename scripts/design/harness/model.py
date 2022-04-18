from dataclasses import dataclass
from torch import nn
import torch


@dataclass
class ModelOutput:
    out: torch.Tensor
    """Raw output of the model."""

    pred: torch.Tensor
    """Prediction of the model. e.g. argmax of the output"""

    loss: torch.Tensor
    """Loss of the model if a loss function was specified."""

    correct: torch.Tensor
    """Correct predictions if labels were specified"""


class Model(nn.Module):
    def model_output(
        self, *, X: torch.Tensor, y: torch.Tensor | None, loss_fn: nn.Module | None
    ) -> ModelOutput:
        """Default implementation for simple models that take the input directly
        and output something the loss function can use directly.

        Args:
            X: Input to the model.
            y: Labels for the input.
            loss_fn: Loss function to use.
        """
        out = self(X)
        pred = out.argmax(dim=1)
        loss = torch.Tensor()
        correct = torch.Tensor()

        if loss_fn:
            if y is None:
                raise ValueError("Computing loss requires labels")
            loss = loss_fn(out, y)

        if y is not None:
            correct = (pred == y).type(torch.float)

        return ModelOutput(out, pred, loss, correct)
