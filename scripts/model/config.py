# These classes mirror the Cfg structs used in memerna.


from dataclasses import dataclass
from enum import Enum
from typing import Optional


class CtdCfg(str, Enum):
    NONE = "none"  #  Do not use CTDs in efn, folding, subopt, partition, etc.
    NO_COAX = "no-coax"  #  Use only terminal mismatches and dangling ends in folding, subopt, partition, etc.
    ALL = "all"  #  Use CTDs in folding, subopt, partition, etc.


@dataclass
class EnergyCfg:
    lonely_pairs: bool = False
    ctd: CtdCfg = CtdCfg.ALL


@dataclass
class SuboptCfg:
    delta: Optional[int] = None
    strucs: Optional[int] = None
    sorted: bool = True