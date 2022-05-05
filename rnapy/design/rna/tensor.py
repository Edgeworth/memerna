import torch


class RnaTensor:
    @staticmethod
    def from_primary(
        primary: str, nuc_map: dict[str, int] = {"X": 0, "A": 1, "C": 2, "G": 3, "U": 4}
    ) -> torch.Tensor:
        """Index representation of a primary structure."""
        return torch.LongTensor([nuc_map[i] for i in primary])

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
