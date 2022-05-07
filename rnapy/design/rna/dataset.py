import random
import RNA
from rnapy.design.rna.tensor import UNK_TOK, RnaTensor
from rnapy.model.random import RandomRna
import torch
from torch.utils.data import Dataset


class RnaDataset(Dataset):
    """Dataset that maps dot-bracket structures to primary structures (designs)"""

    primary_in: list[torch.Tensor]
    primary_out: list[torch.Tensor]
    db: list[torch.Tensor]

    def __init__(self, *, num_struc: int, struc_len: int, max_seq_len: int) -> None:
        if max_seq_len > struc_len:
            raise ValueError(
                f"Transformer input len ({max_seq_len}) > structure size ({struc_len})",
            )

        self.primary_in = []
        self.primary_out = []
        self.db = []
        for _ in range(num_struc):
            primary_out = RandomRna.primary(struc_len)
            primary_in = self.masked_primary(primary_out)
            db, _ = RNA.fold(primary_out)

            out_tensor = RnaTensor.from_primary(primary_out)
            in_tensor = RnaTensor.from_primary(primary_in)
            db_tensor = RnaTensor.from_db(db)
            for i in range(struc_len - max_seq_len):
                self.primary_in.append(in_tensor[i : i + max_seq_len])
                self.primary_out.append(out_tensor[i : i + max_seq_len])
                self.db.append(db_tensor[i : i + max_seq_len])

    def masked_primary(self, primary: str) -> str:
        """Mask off a proportion of the primary sequence"""
        # For now, randomly mask off half.
        indices = random.sample(range(len(primary)), len(primary) // 2)
        bases = list(primary)
        for i in indices:
            bases[i] = UNK_TOK
        return "".join(bases)

    def __len__(self) -> int:
        return len(self.db)

    def __getitem__(self, idx: int) -> list[torch.Tensor]:
        """
        Returns:
            Next value in the sequence for each point in the sequence"""
        return [self.db[idx], self.primary_in[idx], self.primary_out[idx]]
