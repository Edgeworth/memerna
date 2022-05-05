import torch
from torch import nn


class Model(nn.Module):
    def model_prediction(self, out: torch.Tensor) -> torch.Tensor:
        """Make a prediction from the output of the model. e.g. argmax of the output

        Args:
            out: Output of the model
        """
        return out.argmax(dim=1)

    def model_loss(
        self, *, batch: list[torch.Tensor], out: torch.Tensor, loss_fn: nn.Module
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Returns the loss for a given batch and output.

        Args:
            batch: Input batch
            out: Output of the model
            loss_fn: Loss function to use.
        Returns:
            loss: Loss for the model.
            accuracy: Whether predictions were correct or not.
        """
        y = batch[-1]
        return loss_fn(out, y), (self.model_prediction(out) == y).type(torch.float)
