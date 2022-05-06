import torch
from bidict import bidict

NUC_MAP = bidict({"X": 0, "A": 1, "C": 2, "G": 3, "U": 4})


class RnaTensor:
    @staticmethod
    def from_primary(primary: str) -> torch.Tensor:
        """Index representation of a primary structure."""
        return torch.LongTensor([NUC_MAP[i] for i in primary])

    @staticmethod
    def to_primary(index: torch.LongTensor) -> str:
        """Primary structure from index representation."""
        return "".join(
            [NUC_MAP.inverse[i] for i in index]  # pylint: disable=unsubscriptable-object
        )

    @staticmethod
    def secondary_flat_mapping(s: list[int]) -> list[int]:
        """0=>(, 1=>., 2=>)"""
        mapping = []
        for i, p in enumerate(s):
            if p == -1:
                mapping.append(1)
            elif p > i:
                mapping.append(0)
            else:
                mapping.append(2)
        return mapping

    @staticmethod
    def from_secondary(secondary: list[int]) -> torch.Tensor:
        """Index representation of a secondary structure."""
        # TODO: try relative, absolute, and 0/1/2 representations.
        return torch.LongTensor(RnaTensor.secondary_flat_mapping(secondary))
