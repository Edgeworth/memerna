import RNA
from scripts.model.random import RandomRna
import torch
from torch.utils.data import Dataset


class RnaDataset(Dataset):
    """Dataset that maps secondary structures to primary structures (designs)"""

    primaries: list[torch.Tensor]
    secondaries: list[torch.Tensor]

    def __init__(self, *, num_seq: int, seq_size: int) -> None:
        for _ in range(num_seq):
            primary = RandomRna.primary(seq_size)
            db, _ = RNA.fold(primary)
            self.primaries.append(primary)

    def __len__(self) -> int:
        return len(self.primaries)

    def __getitem__(self, idx: int) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Returns:
            Next value in the sequence for each point in the sequence"""
        return (self.secondaries[idx], self.primaries[idx])
