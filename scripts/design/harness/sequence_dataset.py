import torch
from torch.utils.data import Dataset


class SequenceDataset(Dataset):
    """Dataset that produces overlapped sequences from the given data and predicts the next value"""

    seq: torch.Tensor
    seq_size: int

    def __init__(self, *, seq: torch.Tensor, seq_size: int) -> None:
        if not isinstance(seq, torch.Tensor):
            raise TypeError("seq must be a torch.Tensor")
        self.seq = seq
        self.seq_size = seq_size

    def __len__(self) -> int:
        return max(len(self.seq) - self.seq_size, 0)

    def __getitem__(self, idx: int) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Returns:
            Next value in the sequence for each point in the sequence"""
        return (self.seq[idx : idx + self.seq_size], self.seq[idx + 1 : idx + self.seq_size + 1])
