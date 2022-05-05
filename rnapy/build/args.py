# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import cloup
from rnapy.build.config import BuildCfg
from rnapy.build.config import BuildKind
from rnapy.build.config import Compiler
from rnapy.build.config import Sanitizer

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
    cloup.option("--rnastructure/--no-rnastructure", default=False),
    cloup.option("--mpfr/--no-mpfr", default=False),
    cloup.option("--float-bits", type=int, default=64),
)


def build_cfg_from_args(  # pylint: disable=too-many-arguments
    memerna_src_path: Path,
    prefix: Path,
    kind: BuildKind,
    compiler: Compiler,
    sanitizer: Sanitizer,
    mpfr: bool,
    float_bits: int,
    rnastructure: bool,
    iwyu: bool,
    **_kwargs: Any,
) -> BuildCfg:
    return BuildCfg(
        src=memerna_src_path,
        prefix=prefix,
        kind=kind,
        compiler=compiler,
        sanitizer=sanitizer,
        mpfr=mpfr,
        float_bits=float_bits,
        rnastructure=rnastructure,
        iwyu=iwyu,
    )
