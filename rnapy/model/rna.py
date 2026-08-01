# Copyright 2016 Eliot Courtney.
from dataclasses import dataclass
from decimal import Decimal
from typing import override

from rnapy.model.parse.sequence import secondary_to_db, seq_to_primary, validate_secondary


@dataclass
class Rna:
    name: str | None = None
    r: str | None = None
    s: list[int] | None = None
    energy: Decimal | None = None

    def __post_init__(self) -> None:
        if self.r is not None:
            if not self.r:
                raise ValueError("RNA sequence must be non-empty")
            self.r = seq_to_primary(self.r)
        if self.s is not None:
            if not self.s:
                raise ValueError("RNA secondary structure must be non-empty")
            validate_secondary(self.s)
        if self.r is not None and self.s is not None and len(self.r) != len(self.s):
            raise ValueError("primary and secondary must be the same length")
        if self.energy is not None and not self.energy.is_finite():
            raise ValueError("RNA energy must be finite")

    @override
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

    def require_name(self) -> str:
        if not self.name:
            raise ValueError("RNA name must be non-empty")
        return self.name

    def require_sequence(self) -> str:
        if self.r is None:
            raise ValueError(f"RNA {self.name} has no sequence")
        return self.r

    def require_secondary(self) -> list[int]:
        if self.s is None:
            raise ValueError(f"RNA {self.name} has no secondary structure")
        return self.s

    def db(self) -> str:
        return secondary_to_db(self.require_secondary())
