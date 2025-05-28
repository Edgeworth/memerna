# Copyright 2016 Eliot Courtney.
from dataclasses import dataclass
from decimal import Decimal

from rnapy.model.parse.sequence import secondary_to_db


@dataclass
class Rna:
    name: str | None = None
    r: str | None = None
    s: list[int] | None = None
    energy: Decimal | None = None

    def __post_init__(self) -> None:
        if self.r and self.s and len(self.r) != len(self.s):
            raise ValueError("primary and secondary must be the same length")

    def __str__(self) -> str:
        res = ""
        if self.name is not None:
            res += f"{self.name}\n"
        if self.r is not None:
            res += f"{self.r}\n"
        if self.s is not None:
            res += f"{self.db()}\n"
        if self.energy is not None:
            res += f"{self.energy}\n"
        return res

    # Used by command line parsing, for example.
    def __len__(self) -> int:
        if self.r is not None:
            return len(self.r)
        if self.s is not None:
            return len(self.s)
        return 0

    def db(self) -> str:
        if not self.s:
            return ""
        return secondary_to_db(self.s)
