import torch
from torch import nn
import torch.nn.functional as F


class Model(nn.Module):
    def model_prediction(self, out: torch.Tensor) -> torch.Tensor:
        """Make a prediction from the output of the model. e.g. argmax of the output

        Args:
            out: Output of the model
        """
        return out.argmax(dim=1)

    def model_loss(
        self,
        *,
        batch: list[torch.Tensor],
        outs: list[torch.Tensor],
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Returns the loss for a given batch and output.

        Args:
            batch: Input batch
            outs: Outputs of the model
        Returns:
            loss: Loss for the model.
            accuracy: Whether predictions were correct or not.
        """
        y = batch[-1]
        loss = F.cross_entropy(outs[0], y)
        return loss, (self.model_prediction(outs[0]) == y).type(torch.float)
