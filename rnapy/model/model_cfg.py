# Copyright 2022 Eliot Courtney.
# These classes mirror the Cfg structs used in memerna.
from dataclasses import dataclass
from decimal import Decimal
from enum import Enum


class CtdCfg(str, Enum):
    NONE = "none"  #  Do not use CTDs in efn, folding, subopt, partition, etc.
    D2 = "d2"  #  Same as ViennaRNA d2 in efn, folding, subopt, partition, etc.
    #  Use only terminal mismatches and dangling ends in folding, subopt, partition, etc.
    NO_COAX = "no-coax"
    ALL = "all"  #  Use CTDs in folding, subopt, partition, etc.


@dataclass
class EnergyCfg:
    lonely_pairs: bool = False
    ctd: CtdCfg = CtdCfg.ALL


@dataclass
class SuboptCfg:
    delta: Decimal | None = None
    strucs: int | None = None
    sorted: bool = True
