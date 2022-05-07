import RNA
from rnapy.design.rna.config import RnaPipelineConfig
from rnapy.design.rna.tensor import RnaTensor
from rnapy.model.random import RandomRna
import torch
from torch.utils.data import Dataset


class RnaDataset(Dataset):
    """Dataset that maps dot-bracket structures to primary structures (designs)"""

    primary: list[torch.Tensor]
    db: list[torch.Tensor]

    def __init__(self, *, num_struc: int, cfg: RnaPipelineConfig) -> None:
        if cfg.max_seq_len > cfg.struc_len:
            raise ValueError(
                f"Transformer input len ({cfg.max_seq_len}) > structure size ({cfg.struc_len})",
            )

        self.primary = []
        self.db = []
        for _ in range(num_struc):
            primary = RandomRna.primary(cfg.struc_len)
            db, _ = RNA.fold(primary)

            primary_tensor = RnaTensor.from_primary(primary)
            db_tensor = RnaTensor.from_db(db)
            for i in range(cfg.struc_len - cfg.max_seq_len):
                self.primary.append(primary_tensor[i : i + cfg.max_seq_len])
                self.db.append(db_tensor[i : i + cfg.max_seq_len])

    def __len__(self) -> int:
        return len(self.db)

    def __getitem__(self, idx: int) -> list[torch.Tensor]:
        """
        Returns:
            Next value in the sequence for each point in the sequence"""
        return [self.db[idx], self.primary[idx]]
