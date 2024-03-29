import torch
from torch.utils.data import Dataset

from rnapy.design.rna.pipeline_cfg import RnaPipelineCfg
from rnapy.model.random import RandomRna

try:
    from ViennaRNA import RNA
except Exception:
    print("ViennaRNA not installed")


class RnaDataset(Dataset):
    """Dataset that maps dot-bracket structures to primary structures (designs)"""

    primary: list[torch.Tensor]
    db: list[torch.Tensor]

    def __init__(self, *, num_struc: int, cfg: RnaPipelineCfg) -> None:
        self.primary = []
        self.db = []
        for _ in range(num_struc):
            primary, db = self._random_primary_db_pair(cfg.struc_len)

            primary_tensor = cfg.tensor.from_primary(primary)
            db_tensor = cfg.tensor.from_db(db)
            for i in range(len(primary_tensor) - cfg.max_seq_len + 1):
                self.primary.append(primary_tensor[i : i + cfg.max_seq_len])
                self.db.append(db_tensor[i : i + cfg.max_seq_len])

    def _random_primary_db_pair(self, length: int) -> tuple[str, str]:
        primary = RandomRna.primary(length)
        db, _ = RNA.fold(primary)
        return primary, db

    def _overfitting_test_primary_db_pair(self, length: int) -> tuple[str, str]:
        """Map where each base always maps to a specific DB value."""
        primary = RandomRna.primary(length)
        db, _ = RNA.fold(primary)

        overfit_map = {"(": "A", ")": "U", ".": "G"}
        primary = "".join(overfit_map[c] for c in db)
        return primary, db

    def __len__(self) -> int:
        return len(self.db)

    def __getitem__(self, idx: int) -> list[torch.Tensor]:
        """
        Returns:
            Next value in the sequence for each point in the sequence"""
        return [self.db[idx], self.primary[idx]]
