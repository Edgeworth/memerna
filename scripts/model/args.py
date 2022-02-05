# Copyright 2022 Eliot Courtney.
from typing import Optional
from scripts.model.config import CtdCfg, EnergyCfg, SuboptCfg
from scripts.model.parse import db_to_secondary, seq_to_primary
from scripts.model.rna import Rna
import cloup

energy_options = cloup.option_group(
    "Energy options",
    cloup.option(
        "--lonely-pairs/--no-lonely-pairs",
        default=False,
        help="allow lonely pairs",
    ),
    cloup.option(
        "--ctd",
        type=cloup.Choice(list(CtdCfg)),
        default=CtdCfg.ALL,
        help="whether to use CTDs",
    ),
)


def energy_cfg_from_args(lonely_pairs: bool, ctd: CtdCfg, **kwargs):
    return EnergyCfg(lonely_pairs=lonely_pairs, ctd=ctd)


subopt_options = cloup.option_group(
    "Subopt options",
    cloup.option(
        "--subopt-delta",
        type=int,
        help="maximum energy delta from minimum",
    ),
    cloup.option(
        "--subopt-strucs",
        type=int,
        help="maximum number of reported structures",
    ),
    cloup.option(
        "--subopt-sorted/--no-subopt-sorted",
        default=True,
        help="if the structures should be sorted",
    ),
)


def subopt_cfg_from_args(
    subopt_delta: Optional[int], subopt_strucs: Optional[int], subopt_sorted: bool, **kwargs
):
    return SuboptCfg(
        delta=subopt_delta,
        strucs=subopt_strucs,
        sorted=subopt_sorted,
    )
