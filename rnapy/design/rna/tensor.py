import itertools
from typing import Any
from typing import Generator
from typing import Sequence

from bidict import bidict
import torch

BOS_IDX = 0  # beginning of sequence
EOS_IDX = 1  # end of sequence
MASK_IDX = 2  # masked
PAD_IDX = 3  # padding
IDX_END = PAD_IDX  # last special index
PAD_TOK = "_"
MASK_TOK = "X"
NUCS = ["A", "C", "G", "U"]


class RnaTensor:
    def from_primary(self, primary: str) -> torch.Tensor:
        """Index representation of a primary structure."""

    def to_primary(self, index: torch.Tensor | list[int]) -> str:
        """Primary structure from index representation."""

    def primary_dim(self) -> int:
        """Dimension of primary tokens."""

    def from_db(self, db: str) -> torch.Tensor:
        """Index representation of a db structure."""

    def to_db(self, index: torch.Tensor | list[int]) -> str:
        """db structure from index representation."""

    def db_dim(self) -> int:
        """Dimension of db tokens."""


class BasicRnaTensor(RnaTensor):
    db_map: bidict
    primary_map: bidict

    def __init__(self) -> None:
        self.primary_map = bidict(
            {
                ">": BOS_IDX,
                "<": EOS_IDX,
                MASK_TOK: MASK_IDX,
                PAD_TOK: PAD_IDX,
                "A": IDX_END + 1,
                "C": IDX_END + 2,
                "G": IDX_END + 3,
                "U": IDX_END + 4,
            },
        )
        self.db_map = bidict(
            {
                ">": BOS_IDX,
                "<": EOS_IDX,
                MASK_TOK: MASK_IDX,
                PAD_TOK: PAD_IDX,
                "(": IDX_END + 1,
                ".": IDX_END + 2,
                ")": IDX_END + 3,
            },
        )

    def from_primary(self, primary: str) -> torch.Tensor:
        return torch.LongTensor([BOS_IDX] + [self.primary_map[i] for i in primary] + [EOS_IDX])

    def to_primary(self, index: torch.Tensor | list[int]) -> str:
        if isinstance(index, torch.Tensor):
            index = index.tolist()
        return "".join([self.primary_map.inverse[i] for i in index])

    def primary_dim(self) -> int:
        return len(self.primary_map)

    def db_flat_mapping(self, db: str) -> list[int]:
        """Maps (.) directly."""
        return [self.db_map[i] for i in db]

    def from_db(self, db: str) -> torch.Tensor:
        # TODO: try relative, absolute, and 0/1/2 representations.
        return torch.LongTensor([BOS_IDX] + self.db_flat_mapping(db) + [EOS_IDX])

    def to_db(self, index: torch.Tensor | list[int]) -> str:
        if isinstance(index, torch.Tensor):
            index = index.tolist()
        return "".join([self.db_map.inverse[i] for i in index])

    def db_dim(self) -> int:
        return len(self.db_map)


def chunks(seq: Sequence[Any], n: int) -> Generator[Sequence[Any], None, None]:
    for i in range(0, len(seq), n):
        yield seq[i : i + n]


class ChunkedRnaTensor(RnaTensor):
    db_map: bidict
    primary_map: bidict
    chunk_size: int

    def __init__(self, chunk_size: int) -> None:
        self.chunk_size = chunk_size

        self.primary_map = bidict(
            {
                ">": BOS_IDX,
                "<": EOS_IDX,
                MASK_TOK: MASK_IDX,
                PAD_TOK: PAD_IDX,
            },
        )
        self.db_map = bidict(
            {
                ">": BOS_IDX,
                "<": EOS_IDX,
                MASK_TOK: MASK_IDX,
                PAD_TOK: PAD_IDX,
                "(": IDX_END + 1,
                ".": IDX_END + 2,
                ")": IDX_END + 3,
            },
        )
        cur_idx = IDX_END + 1
        iters = []
        for _ in range(chunk_size):
            iters.append(NUCS)
            for chunk in itertools.product(*iters):
                self.primary_map["".join(chunk)] = cur_idx
                cur_idx += 1

    def from_primary(self, primary: str) -> torch.Tensor:
        t = []
        for chunk in chunks(primary, self.chunk_size):
            t.append(self.primary_map["".join(chunk)])

        return torch.LongTensor([BOS_IDX] + t + [EOS_IDX])

    def to_primary(self, index: torch.Tensor | list[int]) -> str:
        if isinstance(index, torch.Tensor):
            index = index.tolist()
        return "".join([self.primary_map.inverse[i] for i in index])

    def primary_dim(self) -> int:
        return len(self.primary_map)

    def from_db(self, db: str) -> torch.Tensor:
        return torch.LongTensor([BOS_IDX] + [self.db_map[i] for i in db] + [EOS_IDX])

    def to_db(self, index: torch.Tensor | list[int]) -> str:
        if isinstance(index, torch.Tensor):
            index = index.tolist()
        return "".join([self.db_map.inverse[i] for i in index])

    def db_dim(self) -> int:
        return len(self.db_map)
