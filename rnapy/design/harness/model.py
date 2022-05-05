import torch
from torch import nn


class Model(nn.Module):
    def model_inputs(self, X: torch.Tensor) -> list[torch.Tensor]:
        """What input the model takes, based on a batch of input."""
        return [X]

    def model_prediction(self, out: torch.Tensor) -> torch.Tensor:
        """Make a prediction from the output of the model. e.g. argmax of the output

        Args:
            out: Output of the model
        """
        return out.argmax(dim=1)

    def model_loss(
        self,
        *,
        out: torch.Tensor,
        y: torch.Tensor,
        loss_fn: nn.Module,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Default loss implementation for simple models that take the input directly
        and output something the loss function can use directly.

        Args:
            out: Output of the model
            y: Labels for the input.
            loss_fn: Loss function to use.
        Returns:
            loss: Loss for the model.
            accuracy: Whether predictions were correct or not.
        """
        return loss_fn(out, y), (self.model_prediction(out) == y).type(torch.float)
