import torch
from bidict import bidict

BOS_IDX = 0  # beginning of sequence
EOS_IDX = 1  # end of sequence
UNK_IDX = 2  # unknown
PRIMARY_MAP = bidict({"X": UNK_IDX, "A": 3, "C": 4, "G": 5, "U": 6})


class RnaTensor:
    @staticmethod
    def from_primary(primary: str) -> torch.Tensor:
        """Index representation of a primary structure."""
        return torch.LongTensor([BOS_IDX] + [PRIMARY_MAP[i] for i in primary] + [EOS_IDX])

    @staticmethod
    def to_primary(index: torch.LongTensor) -> str:
        """Primary structure from index representation."""
        return "".join([PRIMARY_MAP.inverse[i] for i in index if i not in [BOS_IDX, EOS_IDX]])

    @staticmethod
    def primary_dim() -> int:
        """Dimension of primary tokens."""
        return UNK_IDX + 1 + len(PRIMARY_MAP)

    @staticmethod
    def secondary_flat_mapping(s: list[int]) -> list[int]:
        """3=>(, 4=>., 5=>)"""
        mapping = []
        for i, p in enumerate(s):
            if p == -1:
                mapping.append(4)
            elif p > i:
                mapping.append(3)
            else:
                mapping.append(5)
        return mapping

    @staticmethod
    def from_secondary(secondary: list[int]) -> torch.Tensor:
        """Index representation of a secondary structure."""
        # TODO: try relative, absolute, and 0/1/2 representations.
        return torch.LongTensor([BOS_IDX] + RnaTensor.secondary_flat_mapping(secondary) + [EOS_IDX])

    @staticmethod
    def secondary_dim() -> int:
        """Dimension of secondary tokens."""
        return UNK_IDX + 1 + 3
