# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import cloup

from rnapy.build.afl_fuzz import AflFuzzCfg
from rnapy.build.build_cfg import BuildCfg, BuildKind, Compiler, Sanitizer
from rnapy.model.model_cfg import CtdCfg, LonelyPairs
from rnapy.util.util import enum_choice

build_cfg_options = cloup.option_group(
    "Build config options",
    cloup.option(
        "--memerna-src-path",
        type=cloup.Path(exists=True, file_okay=False, resolve_path=True, path_type=Path),
        envvar="MRNA",
        show_envvar=True,
        default=Path(),  # Use the current directory as default.
        help="Path to memerna source directory",
    ),
    cloup.option(
        "--prefix",
        type=cloup.Path(exists=True, file_okay=False, resolve_path=True, path_type=Path),
        default=Path.home() / "bin",
        help="Where to place build directory",
    ),
    cloup.option("--kind", type=enum_choice(BuildKind), default="debug"),
    cloup.option("--compiler", type=enum_choice(Compiler), default="default"),
    cloup.option("--sanitizer", type=enum_choice(Sanitizer), default="none"),
    cloup.option(
        "--iwyu/--no-iwyu", default=False, help="Whether to build with include-what-you-use"
    ),
    cloup.option("--lto/--no-lto", default=False, help="Whether to build with LTO"),
    cloup.option("--rnastructure/--no-rnastructure", default=False),
    cloup.option("--mpfr/--no-mpfr", default=False),
    cloup.option("--float-precision", type=int, default=15),
    cloup.option("--energy-precision", type=int, default=2),
)


def build_cfg_from_args(
    memerna_src_path: Path,
    prefix: Path,
    kind: BuildKind,
    compiler: Compiler,
    sanitizer: Sanitizer,
    mpfr: bool,
    rnastructure: bool,
    iwyu: bool,
    lto: bool,
    float_precision: int,
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
        float_precision=float_precision,
        energy_precision=energy_precision,
    )


afl_fuzz_cfg_options = cloup.option_group(
    "afl-fuzz config options",
    cloup.option("--max-len", type=int, default=None, help="Max sequence length"),
    cloup.option("--random-pf/--no-random-pf", default=None, help="Random pseudofree energies"),
    cloup.option("--energy-model", default=None, help="Energy model"),
    cloup.option("--ctd", type=enum_choice(CtdCfg), default=None, help="CTD mode"),
    cloup.option(
        "--lonely-pairs", type=enum_choice(LonelyPairs), default=None, help="Lonely pairs mode"
    ),
    cloup.option("--backends", multiple=True, help="Backends to fuzz"),
    cloup.option("--brute-max", type=int, default=None, help="Max brute force size"),
    cloup.option("--mfe/--no-mfe", default=None, help="Fuzz MFE"),
    cloup.option(
        "--mfe-rnastructure/--no-mfe-rnastructure", default=None, help="Fuzz MFE RNAstructure"
    ),
    cloup.option("--mfe-table/--no-mfe-table", default=None, help="Check MFE DP tables"),
    cloup.option("--subopt/--no-subopt", default=None, help="Fuzz suboptimal folding"),
    cloup.option(
        "--subopt-rnastructure/--no-subopt-rnastructure",
        default=None,
        help="Fuzz subopt RNAstructure",
    ),
    cloup.option("--subopt-strucs", type=int, default=None, help="Max subopt structures"),
    cloup.option("--subopt-delta", type=float, default=None, help="Max subopt energy delta"),
    cloup.option("--pfn/--no-pfn", default=None, help="Fuzz partition function"),
    cloup.option(
        "--pfn-rnastructure/--no-pfn-rnastructure", default=None, help="Fuzz PFN RNAstructure"
    ),
)

afl_fuzz_index_option = cloup.option(
    "--index",
    type=cloup.IntRange(min=0),
    default=0,
    help="Index of the generated fuzzer configuration to use",
)


def build_afl_fuzz_cfg_from_args(
    build_cfg: BuildCfg,
    max_len: int | None = None,
    random_pf: bool | None = None,
    energy_model: str | None = None,
    ctd: CtdCfg | None = None,
    lonely_pairs: LonelyPairs | None = None,
    backends: tuple[str, ...] = (),
    brute_max: int | None = None,
    mfe: bool | None = None,
    mfe_rnastructure: bool | None = None,
    mfe_table: bool | None = None,
    subopt: bool | None = None,
    subopt_rnastructure: bool | None = None,
    subopt_strucs: int | None = None,
    subopt_delta: float | None = None,
    pfn: bool | None = None,
    pfn_rnastructure: bool | None = None,
    **_kwargs: Any,
) -> AflFuzzCfg:
    return AflFuzzCfg(
        build_cfg=build_cfg,
        fuzz_max_len=max_len,
        fuzz_random_pseudofree=random_pf,
        fuzz_energy_model=energy_model,
        fuzz_ctd=ctd,
        fuzz_lonely_pairs=lonely_pairs,
        fuzz_backends=list(backends) if backends else None,
        fuzz_brute_max=brute_max,
        fuzz_mfe=mfe,
        fuzz_mfe_rnastructure=mfe_rnastructure,
        fuzz_mfe_table=mfe_table,
        fuzz_subopt=subopt,
        fuzz_subopt_rnastructure=subopt_rnastructure,
        fuzz_subopt_strucs=subopt_strucs,
        fuzz_subopt_delta=subopt_delta,
        fuzz_pfn=pfn,
        fuzz_pfn_rnastructure=pfn_rnastructure,
    )
