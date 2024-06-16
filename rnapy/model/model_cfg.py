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

    def desc(self) -> str:
        desc = f"ctd:{self.ctd}-lp:{self.lonely_pairs}"
        if self.energy_model:
            desc += f"-{self.energy_model}"
        if self.backend:
            desc += f"-{self.backend}"
        return desc


@dataclass
class SuboptCfg:
    sorted_strucs: bool
    delta: Decimal | None = None
    strucs: int | None = None
    time_secs: float | None = None
    count_only: bool = False
    # Program specific suboptimal folding algorithm option.
    algorithm: str | None = None

    def desc(self) -> str:
        desc = "subopt"
        if self.delta:
            desc += f"-delta:{self.delta}"
        if self.strucs:
            desc += f"-strucs:{self.strucs}"
        if self.time_secs:
            desc += f"-time:{self.time_secs}"
        if self.sorted_strucs:
            desc += "-sorted"
        if self.algorithm:
            desc += f"-alg:{self.algorithm}"
        return desc
