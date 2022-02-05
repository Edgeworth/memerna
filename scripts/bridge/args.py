# Copyright 2022 Eliot Courtney.
from pathlib import Path
import cloup

from scripts.bridge.memerna import MemeRna
from scripts.bridge.sparsemfefold import SparseMfeFold
from scripts.bridge.unafold import UnaFold


def validate_memerna(ctx, param, value):
    if value is None:
        return None
    return MemeRna(value)


def validate_rnastructure(ctx, param, value):
    if value is None:
        return None
    return RNAstructure(value)


def validate_viennarna(ctx, param, value):
    if value is None:
        return None
    return ViennaRna(value)


def validate_sparsemfefold(ctx, param, value):
    if value is None:
        return None
    return SparseMfeFold(value)


def validate_unafold(ctx, param, value):
    if value is None:
        return None
    return UnaFold(value)


bridge_options = cloup.option_group(
    "Path options",
    cloup.option(
        "--memerna-path",
        "memerna",
        envvar="MRNA",
        show_envvar=True,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
        callback=validate_memerna,
    ),
    cloup.option(
        "--rnastructure-path",
        "rnastructure",
        envvar="RNASTRUCTURE",
        show_envvar=True,
        callback=validate_rnastructure,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
    cloup.option(
        "--viennarna-path",
        "viennarna",
        envvar="VIENNARNA",
        show_envvar=True,
        callback=validate_viennarna,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
    cloup.option(
        "--sparsemfefold-path",
        "sparsemfefold",
        envvar="SPARSEMFEFOLD",
        show_envvar=True,
        callback=validate_sparsemfefold,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
    cloup.option(
        "--unafold-path",
        "unafold",
        envvar="UNAFOLD",
        show_envvar=True,
        callback=validate_unafold,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
)
