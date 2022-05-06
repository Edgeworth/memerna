import RNA
from rnapy.design.rna.tensor import RnaTensor
from rnapy.model.parse.sequence import db_to_secondary
from rnapy.model.random import RandomRna
import torch
from torch.utils.data import Dataset


class RnaDataset(Dataset):
    """Dataset that maps secondary structures to primary structures (designs)"""

    primaries: list[torch.Tensor]
    secondaries: list[torch.Tensor]

    def __init__(self, *, num_struc: int, struc_len: int, max_seq_len: int) -> None:
        if max_seq_len > struc_len:
            raise ValueError(
                f"Transformer input len ({max_seq_len}) > structure size ({struc_len})",
            )

        self.primaries = []
        self.secondaries = []
        for _ in range(num_struc):
            primary = RandomRna.primary(struc_len)
            db, _ = RNA.fold(primary)

            primary_tensor = RnaTensor.from_primary(primary)
            secondary_tensor = RnaTensor.from_secondary(db_to_secondary(db))
            # TODO: providing offset primary here.
            for i in range(struc_len - max_seq_len - 1):
                self.primaries.append(primary_tensor[i + 1 : i + max_seq_len + 1])
                self.secondaries.append(secondary_tensor[i : i + max_seq_len])

    def __len__(self) -> int:
        return len(self.primaries)

    def __getitem__(self, idx: int) -> list[torch.Tensor]:
        """
        Returns:
            Next value in the sequence for each point in the sequence"""
        return [self.secondaries[idx], self.primaries[idx]]
