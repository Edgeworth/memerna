# Copyright 2022 Eliot Courtney.
from enum import StrEnum
from pathlib import Path

import click
import cloup

from rnapy.bridge.linearfold import LinearFold
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.memerna01 import MemeRna01
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.sparsemfefold import SparseMFEFold
from rnapy.bridge.sparsernafold import SparseRNAFolD
from rnapy.bridge.viennarna import ViennaRna


class BridgeName(StrEnum):
    LINEARFOLD = "linearfold"
    MEMERNA = "memerna"
    MEMERNA01 = "memerna01"
    RNASTRUCTURE = "rnastructure"
    SPARSEMFEFOLD = "sparsemfefold"
    SPARSERNAFOLD = "sparsernafold"
    VIENNARNA = "viennarna"


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


memerna_path_option = cloup.option(
    "--memerna-path",
    BridgeName.MEMERNA.value,
    envvar="MRNA_DIST",
    show_envvar=True,
    type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
    callback=validate_memerna,
    help="path to memerna build directory",
)

memerna01_path_option = cloup.option(
    "--memerna01-path",
    BridgeName.MEMERNA01.value,
    envvar="MRNA01_DIST",
    show_envvar=True,
    type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
    callback=validate_memerna01,
    help="path to memerna0.1 build directory",
)

linearfold_path_option = cloup.option(
    "--linearfold-path",
    BridgeName.LINEARFOLD.value,
    envvar="LINEARFOLD",
    show_envvar=True,
    callback=validate_linearfold,
    type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
)

rnastructure_path_option = cloup.option(
    "--rnastructure-path",
    BridgeName.RNASTRUCTURE.value,
    envvar="RNASTRUCTURE",
    show_envvar=True,
    callback=validate_rnastructure,
    type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
    help="path to RNAstructure source directory (in-tree build assumed to be done)",
)

sparsemfefold_path_option = cloup.option(
    "--sparsemfefold-path",
    BridgeName.SPARSEMFEFOLD.value,
    envvar="SPARSEMFEFOLD",
    show_envvar=True,
    callback=validate_sparsemfefold,
    type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
)

viennarna_path_option = cloup.option(
    "--viennarna-path",
    BridgeName.VIENNARNA.value,
    envvar="VIENNARNA",
    show_envvar=True,
    callback=validate_viennarna,
    type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
)

sparsernafold_path_option = cloup.option(
    "--sparsernafold-path",
    BridgeName.SPARSERNAFOLD.value,
    envvar="SPARSERNAFOLD",
    show_envvar=True,
    callback=validate_sparsernafold,
    type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
)

bridge_options = cloup.option_group(
    "Bridge options",
    memerna_path_option,
    memerna01_path_option,
    linearfold_path_option,
    rnastructure_path_option,
    sparsemfefold_path_option,
    viennarna_path_option,
    sparsernafold_path_option,
)
