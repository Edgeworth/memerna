# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import cloup
from rnapy.build.afl_fuzz import AflFuzzCfg
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

afl_fuzz_cfg_options = cloup.option_group(
    "afl-fuzz config options",
    cloup.option(
        "--max-len",
        default=-1,
        help="Max size of sequences to fuzz",
    ),
    cloup.option(
        "--random-model/--no-random-model",
        default=False,
        help="Whether to use random model for energy",
    ),
    cloup.option(
        "--brute-max",
        default=22,
        help="Max size of sequences to brute force",
    ),
    cloup.option(
        "--mfe/--no-mfe",
        default=False,
        help="Whether to fuzz mfe",
    ),
    cloup.option(
        "--mfe-rnastructure/--no-mfe-rnastructure",
        default=False,
        help="Whether to fuzz mfe with rnastructure",
    ),
    cloup.option(
        "--mfe-table/--no-mfe-table",
        default=False,
        help="Whether to fuzz mfe with table",
    ),
    cloup.option(
        "--subopt/--no-subopt",
        default=False,
        help="Whether to fuzz subopt",
    ),
    cloup.option(
        "--subopt-rnastructure/--no-subopt-rnastructure",
        default=False,
        help="Whether to fuzz subopt with rnastructure",
    ),
    cloup.option(
        "--subopt-strucs",
        default=5000,
        help="Maximum number of structures to generate for subopt",
    ),
    cloup.option(
        "--subopt-delta",
        default=0.6,
        help="Maximum energy delta for subopt",
    ),
    cloup.option(
        "--part/--no-part",
        default=False,
        help="Whether to fuzz partition function",
    ),
    cloup.option(
        "--part-rnastructure/--no-part-rnastructure",
        default=False,
        help="Whether to fuzz partition function with rnastructure",
    ),
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


def build_afl_fuzz_cfg_from_args(  # pylint: disable=too-many-arguments
    build_cfg: BuildCfg,
    max_len: int = -1,
    random_model: bool = False,
    brute_max: int = 22,
    mfe: bool = False,
    mfe_rnastructure: bool = False,
    mfe_table: bool = False,
    subopt: bool = False,
    subopt_rnastructure: bool = False,
    subopt_strucs: int = 5000,
    subopt_delta: float = 0.6,
    part: bool = False,
    part_rnastructure: bool = False,
    **_kwargs: Any,
) -> AflFuzzCfg:
    return AflFuzzCfg(
        build_cfg=build_cfg,
        fuzz_max_len=max_len,
        fuzz_random_model=random_model,
        fuzz_brute_max=brute_max,
        fuzz_mfe=mfe,
        fuzz_mfe_rnastructure=mfe_rnastructure,
        fuzz_mfe_table=mfe_table,
        fuzz_subopt=subopt,
        fuzz_subopt_rnastructure=subopt_rnastructure,
        fuzz_subopt_strucs=subopt_strucs,
        fuzz_subopt_delta=subopt_delta,
        fuzz_part=part,
        fuzz_part_rnastructure=part_rnastructure,
    )
