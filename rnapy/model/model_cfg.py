# Copyright 2022 Eliot Courtney.
# These classes mirror the Cfg structs used in memerna.
from dataclasses import dataclass
from decimal import Decimal
from enum import StrEnum


class LonelyPairs(StrEnum):
    OFF = "off"
    HEURISTIC = "heuristic"
    ON = "on"


class CtdCfg(StrEnum):
    NONE = "none"  #  Do not use CTDs in efn, folding, subopt, partition, etc.
    D2 = "d2"  #  Same as ViennaRNA d2 in efn, folding, subopt, partition, etc.
    #  Use only terminal mismatches and dangling ends in folding, subopt, partition, etc.
    NO_COAX = "no-coax"
    ALL = "all"  #  Use CTDs in folding, subopt, partition, etc.


@dataclass
class EnergyCfg:
    ctd: CtdCfg
    lonely_pairs: LonelyPairs
    # Program specific energy model option.
    energy_model: str | None = None
    backend: str | None = None

    def desc(self) -> dict[str, str]:
        return {
            "ctd": self.ctd.value,
            "lonely_pairs": self.lonely_pairs.value,
            "energy_model": self.energy_model or "",
            "backend": self.backend or "",
        }


@dataclass
class SuboptCfg:
    sorted_strucs: bool
    delta: Decimal | None = None
    strucs: int | None = None
    time_secs: float | None = None
    count_only: bool = False
    # Program specific suboptimal folding algorithm option.
    algorithm: str | None = None

    def desc(self) -> dict[str, str]:
        return {
            "sorted_strucs": str(self.sorted_strucs),
            "delta": str(self.delta) if self.delta is not None else "",
            "strucs": str(self.strucs) if self.strucs is not None else "",
            "time_secs": str(self.time_secs) if self.time_secs is not None else "",
            "count_only": str(self.count_only),
            "algorithm": self.algorithm or "",
        }
