# Copyright 2022 E.
from pathlib import Path
from typing import Any

import cloup
from rnapy.build.build_cfg import BuildCfg
from rnapy.build.build_cfg import BuildKind
from rnapy.build.build_cfg import Compiler
from rnapy.build.build_cfg import Sanitizer

build_cfg_options = cloup.option_group(
    "Build config options",
    cloup.option(
        "--memerna-src-path",
        type=cloup.Path(exists=True, file_okay=False, resolve_path=True, path_type=Path),
        envvar="MRNA",
        show_envvar=True,
        required=True,
        help="Path to memerna source directory",
    ),
    cloup.option(
        "--prefix",
        type=cloup.Path(exists=True, file_okay=False, resolve_path=True, path_type=Path),
        default=Path.home() / "bin",
        help="Where to place build directory",
    ),
    cloup.option(
        "--kind",
        type=cloup.Choice(list(BuildKind)),
        default="debug",
    ),
    cloup.option(
        "--compiler",
        type=cloup.Choice(list(Compiler)),
        default="default",
    ),
    cloup.option(
        "--sanitizer",
        type=cloup.Choice(list(Sanitizer)),
        default="none",
    ),
    cloup.option(
        "--iwyu/--no-iwyu",
        default=False,
        help="Whether to build with include-what-you-use",
    ),
    cloup.option(
        "--lto/--no-lto",
        default=False,
        help="Whether to build with LTO",
    ),
    cloup.option("--rnastructure/--no-rnastructure", default=False),
    cloup.option("--mpfr/--no-mpfr", default=False),
    cloup.option("--float-bits", type=int, default=64),
    cloup.option("--energy-precision", type=int, default=2),
)


def build_cfg_from_args(  # pylint: disable=too-many-arguments
    memerna_src_path: Path,
    prefix: Path,
    kind: BuildKind,
    compiler: Compiler,
    sanitizer: Sanitizer,
    mpfr: bool,
    rnastructure: bool,
    iwyu: bool,
    lto: bool,
    float_bits: int,
    energy_precision: int,
    **_kwargs: Any,
) -> BuildCfg:
    return BuildCfg(
        src=memerna_src_path,
        prefix=prefix,
        kind=kind,
        compiler=compiler,
        sanitizer=sanitizer,
        mpfr=mpfr,
        rnastructure=rnastructure,
        iwyu=iwyu,
        lto=lto,
        float_bits=float_bits,
        energy_precision=energy_precision,
    )
