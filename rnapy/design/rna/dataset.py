import RNA
from rnapy.design.rna.tensor import RnaTensor
from rnapy.model.random import RandomRna
import torch
from torch.utils.data import Dataset


class RnaDataset(Dataset):
    """Dataset that maps dot-bracket structures to primary structures (designs)"""

    primaries: list[torch.Tensor]
    dbs: list[torch.Tensor]

    def __init__(self, *, num_struc: int, struc_len: int, max_seq_len: int) -> None:
        if max_seq_len > struc_len:
            raise ValueError(
                f"Transformer input len ({max_seq_len}) > structure size ({struc_len})",
            )

        self.primaries = []
        self.dbs = []
        for _ in range(num_struc):
            primary = RandomRna.primary(struc_len)
            db, _ = RNA.fold(primary)

            primary_tensor = RnaTensor.from_primary(primary)
            db_tensor = RnaTensor.from_db(db)
            # TODO: Includes an extra nuc from the primary sequence here.
            # For now, the transformer predicts the next token, so the
            # input to the decodeer is [i:i+max_seq_len] and the output is
            # trained to be [i+1:i+max_seq_len+1]
            for i in range(struc_len - max_seq_len - 1):
                self.primaries.append(primary_tensor[i : i + max_seq_len + 1])
                self.dbs.append(db_tensor[i : i + max_seq_len])

    def __len__(self) -> int:
        return len(self.primaries)

    def __getitem__(self, idx: int) -> list[torch.Tensor]:
        """
        Returns:
            Next value in the sequence for each point in the sequence"""
        return [self.dbs[idx], self.primaries[idx]]
