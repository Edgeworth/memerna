# Copyright 2022 Eliot Courtney.

import logging
from typing import Any

import click
import cloup

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.util.command import CmdLimits


def validate_cpu_affinity(
    ctx: click.Context, param: click.Parameter, value: str | None
) -> set[int] | None:
    if value is None:
        return None

    cpu_affinity: set[int] = set()
    parts = [part.strip() for part in value.strip().split(",")]
    for part in parts:
        try:
            cpu_id = int(part)
            if cpu_id < 0:
                raise click.BadParameter(
                    f"CPU ID '{part}' must be non-negative.", ctx=ctx, param=param
                )
            cpu_affinity.add(cpu_id)
        except ValueError as ex:
            raise click.BadParameter(
                f"CPU ID '{part}' must be an integer.", ctx=ctx, param=param
            ) from ex

    return cpu_affinity


limit_options = cloup.option_group(
    "Resource limit options",
    cloup.option("--time-limit-seconds", type=int, help="maximum time to run any bridge packages"),
    cloup.option(
        "--memory-limit-bytes", type=int, help="maximum memory to use for any bridge packages"
    ),
    cloup.option(
        "--cpu-affinity",
        type=str,
        callback=validate_cpu_affinity,
        help="Set CPU affinity (e.g., '0', '0,1,2').",
    ),
    help="Constraints for command execution",
)


def cmd_limits_from_args(
    time_limit_seconds: int | None,
    memory_limit_bytes: int | None,
    cpu_affinity: set[int] | None,
    **_kwargs: Any,
) -> CmdLimits:
    return CmdLimits(
        time_sec=time_limit_seconds, mem_bytes=memory_limit_bytes, cpu_affinity=cpu_affinity
    )


def init_package_limits(**kwargs: Any) -> None:
    limits = cmd_limits_from_args(**kwargs)
    for v in kwargs.values():
        if isinstance(v, RnaPackage):
            logging.getLogger(__name__).info(
                f"Setting limits on package {v.package_name()}: {limits}"
            )
            v.limits = limits
