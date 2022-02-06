#!/usr/bin/env python3
# Copyright 2016 E.
import os
from pathlib import Path
import shutil
import sys

import click
import cloup


def run_command(cmd):
    if os.system(cmd) != 0:
        sys.exit(1)


def get_compilers(compiler):
    if compiler == "default":
        return ("cc", "c++")
    elif compiler == "gcc":
        return ("gcc", "g++")
    elif compiler == "clang":
        return ("clang", "clang++")
    elif compiler == "afl-lto":
        return ("afl-clang-lto", "afl-clang-lto++")
    elif compiler == "afl-fast":
        return ("afl-clang-fast", "afl-clang-fast++")
    else:
        raise f"Unknown compiler {compiler}"


def get_sanitizers(sanitizer):
    if sanitizer == "none":
        return {}
    elif sanitizer == "asan":
        return {"SANITIZE_ADDRESS": "ON"}
    elif sanitizer == "tsan":
        return {"SANITIZE_THREAD": "ON"}
    elif sanitizer == "ubsan":
        return {"SANITIZE_UNDEFINED": "ON"}
    else:
        raise f"Unknown sanitizer {sanitizer}"


@cloup.command()
# Build environment options
@cloup.option(
    "--prefix",
    type=Path,
    default=os.path.join(Path.home(), "bin"),
    help="Where to place build directory",
)
@cloup.option(
    "--memerna-src-path",
    type=cloup.Path(exists=True, resolve_path=True, file_okay=False, path_type=Path),
    envvar="MRNA",
    show_envvar=True,
    required=True,
    help="Path to memerna source directory",
)
@cloup.option(
    "--type",
    type=cloup.Choice(["debug", "release", "relwithdebinfo"]),
    default="debug",
)
@cloup.option(
    "--compiler",
    type=cloup.Choice(["default", "gcc", "clang", "afl-fast", "afl-lto"]),
    default="default",
)
@cloup.option(
    "--sanitizer",
    type=cloup.Choice(["none", "asan", "tsan", "ubsan"]),
    default="none",
)
@cloup.option("--iwyu/--no-iwyu", default=False, help="Whether to build with include-what-you-use")
@cloup.option("--regenerate/--no-regenerate", default=False)
@cloup.argument("targets", nargs=-1)
# Memerna configuration options:
@cloup.option("--rnastructure/--no-rnastructure", default=False)
@cloup.option("--mpfr/--no-mpfr", default=False)
@cloup.option("--float-bits", type=int, default=64)
# Misc options:
@cloup.option("--test/--no-test", default=False)
@cloup.option("--build/--no-build", default=True)
def build(
    prefix,
    memerna_src_path,
    type,
    compiler,
    sanitizer,
    iwyu,
    regenerate,
    targets,
    rnastructure,
    mpfr,
    float_bits,
    test,
    build,
):
    env = [
        # Add stack protector etc to catch non-crashing memory bugs.
        "AFL_HARDEN=1",
    ]
    compilers = get_compilers(compiler)

    defs = {
        "CMAKE_C_COMPILER": compilers[0],
        "CMAKE_CXX_COMPILER": compilers[1],
        "CMAKE_BUILD_TYPE": type,
        "USE_RNASTRUCTURE": "ON" if rnastructure else "OFF",
        "USE_MPFR": "ON" if mpfr else "OFF",
        "USE_IWYU": "ON" if iwyu else "OFF",
        "FLOAT_BITS": f"{float_bits}",
    }
    defs.update(get_sanitizers(sanitizer))

    build_dir = os.path.join(
        prefix,
        "memerna",
        defs["CMAKE_CXX_COMPILER"] + "-" + defs["CMAKE_BUILD_TYPE"],
    )
    if rnastructure:
        build_dir += "-rnastructure"
    if mpfr:
        build_dir += "-mpfr"
    build_dir += f"-{float_bits}"
    if iwyu:
        build_dir += "-iwyu"

    if sanitizer != "none":
        build_dir += f"-{sanitizer}"

    if regenerate and os.path.exists(build_dir):
        shutil.rmtree(build_dir)

    click.echo(build_dir)
    if not os.path.exists(build_dir):
        os.makedirs(build_dir)
        regenerate = True

    os.chdir(build_dir)
    if regenerate:
        click.echo("Regenerating cmake files.")
        def_str = " ".join(f"-D {i}={k}" for i, k in defs.items())
        run_command(f"cmake {def_str} {memerna_src_path}")

    if build:
        run_command(f"{' '.join(env)} make -j$(($(nproc)-1)) {' '.join(targets)}")

    if test:
        run_command(f"./run_tests")

    os.chdir(memerna_src_path)


if __name__ == "__main__":
    build()
