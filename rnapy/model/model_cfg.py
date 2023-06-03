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
    lonely_pairs: LonelyPairs = LonelyPairs.HEURISTIC
    ctd: CtdCfg = CtdCfg.ALL
    # Program specific energy model option.
    model: str | None = None


@dataclass
class SuboptCfg:
    delta: Decimal | None = None
    strucs: int | None = None
    sorted: bool = True
