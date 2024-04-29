# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import click
import cloup

from rnapy.bridge.linearfold import LinearFold
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.memerna01 import MemeRna01
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.sparsemfefold import SparseMFEFold
from rnapy.bridge.sparsernafold import SparseRNAFolD
from rnapy.bridge.viennarna import ViennaRna


def validate_memerna(
    _ctx: click.Context, _param: click.Parameter, value: Path | None
) -> MemeRna | None:
    if value is None:
        return None
    return MemeRna(value)


def validate_memerna01(
    _ctx: click.Context, _param: click.Parameter, value: Path | None
) -> MemeRna01 | None:
    if value is None:
        return None
    return MemeRna01(value)


def validate_linearfold(
    _ctx: click.Context, _param: click.Parameter, value: Path | None
) -> LinearFold | None:
    if value is None:
        return None
    return LinearFold(value)


def validate_rnastructure(
    _ctx: click.Context, _param: click.Parameter, value: Path | None
) -> RNAstructure | None:
    if value is None:
        return None
    return RNAstructure(value)


def validate_sparsemfefold(
    _ctx: click.Context, _param: click.Parameter, value: Path | None
) -> SparseMFEFold | None:
    if value is None:
        return None
    return SparseMFEFold(value)


def validate_viennarna(
    _ctx: click.Context, _param: click.Parameter, value: Path | None
) -> ViennaRna | None:
    if value is None:
        return None
    return ViennaRna(value)


def validate_sparsernafold(
    _ctx: click.Context, _param: click.Parameter, value: Path | None
) -> SparseRNAFolD | None:
    if value is None:
        return None
    return SparseRNAFolD(value)


bridge_options = cloup.option_group(
    "Bridge options",
    cloup.option("--time-limit-seconds", type=int, help="maximum time to run any bridge packages"),
    cloup.option(
        "--memory-limit-bytes", type=int, help="maximum memory to use for any bridge packages"
    ),
    cloup.option(
        "--memerna-path",
        "memerna",
        envvar="MRNA_DIST",
        show_envvar=True,
        type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
        callback=validate_memerna,
        help="path to memerna build directory",
    ),
    cloup.option(
        "--memerna01-path",
        "memerna01",
        envvar="MRNA01_DIST",
        show_envvar=True,
        type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
        callback=validate_memerna01,
        help="path to memerna0.1 build directory",
    ),
    cloup.option(
        "--linearfold-path",
        "linearfold",
        envvar="LINEARFOLD",
        show_envvar=True,
        callback=validate_linearfold,
        type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
    ),
    cloup.option(
        "--rnastructure-path",
        "rnastructure",
        envvar="RNASTRUCTURE",
        show_envvar=True,
        callback=validate_rnastructure,
        type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
        help="path to RNAstructure source directory (in-tree build assumed to be done)",
    ),
    cloup.option(
        "--sparsemfefold-path",
        "sparsemfefold",
        envvar="SPARSEMFEFOLD",
        show_envvar=True,
        callback=validate_sparsemfefold,
        type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
        help="",
    ),
    cloup.option(
        "--viennarna-path",
        "viennarna",
        envvar="VIENNARNA",
        show_envvar=True,
        callback=validate_viennarna,
        type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
    ),
    cloup.option(
        "--sparsernafold-path",
        "sparsernafold",
        envvar="SPARSERNAFOLD",
        show_envvar=True,
        callback=validate_sparsernafold,
        type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
    ),
)


def init_package_limits(
    time_limit_seconds: int | None, memory_limit_bytes: int | None, **kwargs: Any
) -> None:
    for v in kwargs.values():
        if isinstance(v, RnaPackage):
            v.limits.time_sec = time_limit_seconds
            v.limits.mem_bytes = memory_limit_bytes
