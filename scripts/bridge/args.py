# Copyright 2022 Eliot Courtney.
from pathlib import Path
import cloup

from scripts.bridge.memerna import MemeRna


def validate_memerna(ctx, param, value):
    if value is None:
        return None
    return MemeRna(value)


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
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
    cloup.option(
        "--viennarna-path",
        "viennarna",
        envvar="VIENNARNA",
        show_envvar=True,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
    cloup.option(
        "--sjsviennarna-path",
        "sjsviennarna",
        envvar="SJSVIENNARNA",
        show_envvar=True,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
    cloup.option(
        "--sjsviennarnampi-path",
        "sjsviennarnampi",
        envvar="SJSVIENNARNAMPI",
        show_envvar=True,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
    cloup.option(
        "--sparsemfefold-path",
        "sparsemfefold",
        envvar="SPARSEMFEFOLD",
        show_envvar=True,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
    cloup.option(
        "--unafold-path",
        "unafold",
        envvar="UNAFOLD",
        show_envvar=True,
        type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    ),
)
