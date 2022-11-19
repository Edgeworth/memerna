# Copyright 2022 Eliot Courtney.
from decimal import Decimal
from typing import Any

from click_params import DECIMAL
import cloup
from rnapy.model.model_cfg import CtdCfg
from rnapy.model.model_cfg import EnergyCfg
from rnapy.model.model_cfg import SuboptCfg

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


def energy_cfg_from_args(lonely_pairs: bool, ctd: CtdCfg, **_kwargs: Any) -> EnergyCfg:
    return EnergyCfg(lonely_pairs=lonely_pairs, ctd=ctd)


subopt_options = cloup.option_group(
    "Subopt options",
    cloup.option(
        "--subopt-delta",
        type=DECIMAL,
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
    subopt_delta: Decimal | None,
    subopt_strucs: int | None,
    subopt_sorted: bool,
    **_kwargs: Any,
) -> SuboptCfg:
    return SuboptCfg(
        delta=subopt_delta,
        strucs=subopt_strucs,
        sorted=subopt_sorted,
    )
