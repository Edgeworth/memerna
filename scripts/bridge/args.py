# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Optional
import cloup

from scripts.bridge.memerna import MemeRna
from scripts.bridge.rnapackage import RnaPackage
from scripts.bridge.sparsemfefold import SparseMfeFold


def validate_memerna(ctx, param, value):
    if value is None:
        return None
    return MemeRna(value)


def validate_rnastructure(ctx, param, value):
    if value is None:
        return None
    return RNAstructure(value)


def validate_sparsemfefold(ctx, param, value):
    if value is None:
        return None
    return SparseMfeFold(value)


def validate_viennarna(ctx, param, value):
    if value is None:
        return None
    return ViennaRna(value)


bridge_options = cloup.option_group(
    "Bridge options",
    cloup.option(
        "--time-limit-seconds",
        type=int,
        help="maximum time to run any bridge packages",
    ),
    cloup.option(
        "--memory-limit-bytes",
        type=int,
        help="maximum memory to use for any bridge packages",
    ),
    cloup.option(
        "--memerna-path",
        "memerna",
        envvar="MRNA",
        show_envvar=True,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
        callback=validate_memerna,
        help="path to memerna build directory"
    ),
    cloup.option(
        "--rnastructure-path",
        "rnastructure",
        envvar="RNASTRUCTURE",
        show_envvar=True,
        callback=validate_rnastructure,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
        help="path to RNAstructure source directory (in-tree build assumed to be done)"
    ),
    cloup.option(
        "--sparsemfefold-path",
        "sparsemfefold",
        envvar="SPARSEMFEFOLD",
        show_envvar=True,
        callback=validate_sparsemfefold,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
        help=""
    ),
    cloup.option(
        "--viennarna-path",
        "viennarna",
        envvar="VIENNARNA",
        show_envvar=True,
        callback=validate_viennarna,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
)


def init_package_limits(
    time_limit_seconds: Optional[int], memory_limit_bytes: Optional[int], **kwargs
):
    for v in kwargs.values():
        if isinstance(v, RnaPackage):
            v.limits.time_sec = time_limit_seconds
            v.limits.mem_bytes = memory_limit_bytes
