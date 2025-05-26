# Copyright 2022 Eliot Courtney.
from decimal import Decimal
from typing import Any

import cloup
from click_params import DECIMAL

from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs, SuboptCfg
from rnapy.util.util import enum_choice

energy_options = cloup.option_group(
    "Energy options",
    cloup.option(
        "--lonely-pairs",
        type=enum_choice(LonelyPairs),
        default=LonelyPairs.HEURISTIC,
        help="allow lonely pairs",
    ),
    cloup.option("--ctd", type=enum_choice(CtdCfg), default=CtdCfg.ALL, help="whether to use CTDs"),
    cloup.option(
        "--energy-model",
        type=str,
        required=False,
        help="energy model to use (only applies to some programs, like memerna)",
    ),
    cloup.option(
        "--backend",
        type=str,
        required=False,
        help="backend to use (only applies to some programs, like memerna)",
    ),
)


def energy_cfg_from_args(
    lonely_pairs: LonelyPairs,
    ctd: CtdCfg,
    energy_model: str | None,
    backend: str | None,
    **_kwargs: Any,
) -> EnergyCfg:
    return EnergyCfg(lonely_pairs=lonely_pairs, ctd=ctd, energy_model=energy_model, backend=backend)


subopt_options = cloup.option_group(
    "Subopt options",
    cloup.option("--subopt-delta", type=DECIMAL, help="maximum energy delta from minimum"),
    cloup.option("--subopt-strucs", type=int, help="maximum number of reported structures"),
    cloup.option(
        "--subopt-time-secs",
        type=float,
        help="maximum time in seconds to run suboptimal folding for",
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
    subopt_time_secs: float | None,
    subopt_sorted: bool,
    **_kwargs: Any,
) -> SuboptCfg:
    return SuboptCfg(
        delta=subopt_delta,
        strucs=subopt_strucs,
        time_secs=subopt_time_secs,
        sorted_strucs=subopt_sorted,
    )
