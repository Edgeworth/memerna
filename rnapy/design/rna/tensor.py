import torch
from bidict import bidict

BOS_IDX = 0  # beginning of sequence
EOS_IDX = 1  # end of sequence
MASK_IDX = 2  # masked
PAD_IDX = 3  # padding
IDX_END = PAD_IDX  # last special index
MASK_TOK = "X"
PRIMARY_MAP = bidict({MASK_TOK: MASK_IDX, "A": 3, "C": 4, "G": 5, "U": 6})
DB_MAP = bidict({MASK_TOK: MASK_IDX, "(": 3, ".": 4, ")": 5})


class RnaTensor:
    @staticmethod
    def from_primary(primary: str) -> torch.Tensor:
        """Index representation of a primary structure."""
        return torch.LongTensor([BOS_IDX] + [PRIMARY_MAP[i] for i in primary] + [EOS_IDX])

    @staticmethod
    def to_primary(index: torch.Tensor | list[int]) -> str:
        """Primary structure from index representation."""
        if isinstance(index, torch.Tensor):
            index = index.tolist()
        return "".join([PRIMARY_MAP.inverse[i] for i in index if i not in [BOS_IDX, EOS_IDX]])

    @staticmethod
    def primary_dim() -> int:
        """Dimension of primary tokens."""
        return IDX_END + 1 + len(PRIMARY_MAP)

    @staticmethod
    def db_flat_mapping(db: str) -> list[int]:
        """3=>(, 4=>., 5=>)"""
        return [DB_MAP[i] for i in db]

    @staticmethod
    def from_db(db: str) -> torch.Tensor:
        """Index representation of a db structure."""
        # TODO: try relative, absolute, and 0/1/2 representations.
        return torch.LongTensor([BOS_IDX] + RnaTensor.db_flat_mapping(db) + [EOS_IDX])

    @staticmethod
    def to_db(index: torch.Tensor | list[int]) -> str:
        """db structure from index representation."""
        if isinstance(index, torch.Tensor):
            index = index.tolist()
        return "".join([DB_MAP.inverse[i] for i in index if i not in [BOS_IDX, EOS_IDX]])

    @staticmethod
    def db_dim() -> int:
        """Dimension of db tokens."""
        return IDX_END + 1 + len(DB_MAP)
