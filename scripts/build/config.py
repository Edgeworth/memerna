# Copyright 2022 E.
from dataclasses import dataclass
from dataclasses import field
from enum import Enum
from pathlib import Path
import shutil

import click
from scripts.util.command import run_shell


class Compiler(str, Enum):
    DEFAULT = "default"
    GCC = "gcc"
    CLANG = "clang"
    AFL_LTO = "afl-lto"
    AFL_FAST = "afl-fast"

    def cc_cxx(self) -> tuple[str, str]:
        if self == Compiler.DEFAULT:
            return "cc", "c++"
        if self == Compiler.GCC:
            return "gcc", "g++"
        if self == Compiler.CLANG:
            return "clang", "clang++"
        if self == Compiler.AFL_LTO:
            return "afl-clang-lto", "afl-clang-lto++"
        if self == Compiler.AFL_FAST:
            return "afl-clang-fast", "afl-clang-fast++"
        raise ValueError(f"Unknown compiler {self}")


class Sanitizer(str, Enum):
    NONE = "none"
    ASAN = "asan"
    TSAN = "tsan"
    UBSAN = "ubsan"

    def cmake_defs(self) -> dict[str, str]:
        if self == Sanitizer.NONE:
            return {}
        if self == Sanitizer.ASAN:
            return {"SANITIZE_ADDRESS": "ON"}
        if self == Sanitizer.TSAN:
            return {"SANITIZE_THREAD": "ON"}
        if self == Sanitizer.UBSAN:
            return {"SANITIZE_UNDEFINED": "ON"}
        raise ValueError(f"Unknown sanitizer {self}")


class BuildKind(str, Enum):
    DEBUG = "debug"
    RELEASE = "release"
    RELWITHDEBINFO = "relwithdebinfo"


@dataclass
class BuildCfg:
    src: Path
    prefix: Path
    env: dict[str, str] = field(default_factory=dict)
    kind: BuildKind = BuildKind.DEBUG
    compiler: Compiler = Compiler.DEFAULT
    sanitizer: Sanitizer = Sanitizer.NONE
    mpfr: bool = False
    float_bits: int = 64
    rnastructure: bool = False
    iwyu: bool = False

    def ident(self) -> str:
        ident = self.kind + "-" + self.compiler
        if self.sanitizer != Sanitizer.NONE:
            ident += f"-{self.sanitizer}"
        if self.mpfr:
            ident += "-mpfr"
        ident += f"-{self.float_bits}"
        if self.rnastructure:
            ident += "-rnastructure"
        if self.iwyu:
            ident += "-iwyu"
        if self.env:
            ident += "-" + "-".join(f"{k}-{v}" for k, v in self.env.items())
        return ident

    def build_path(self) -> Path:
        return self.prefix / "memerna" / self.ident()

    def is_afl(self) -> bool:
        return self.compiler in [Compiler.AFL_LTO, Compiler.AFL_FAST]

    def _generate_cmake(self) -> None:
        cc, cxx = self.compiler.cc_cxx()

        defs = {
            "CMAKE_C_COMPILER": cc,
            "CMAKE_CXX_COMPILER": cxx,
            "CMAKE_BUILD_TYPE": self.kind,
            "USE_RNASTRUCTURE": "ON" if self.rnastructure else "OFF",
            "USE_MPFR": "ON" if self.mpfr else "OFF",
            "USE_IWYU": "ON" if self.iwyu else "OFF",
            "FLOAT_BITS": f"{self.float_bits}",
        }
        defs.update(self.sanitizer.cmake_defs())
        def_str = " ".join(f"-D {i}={k}" for i, k in defs.items())

        build_path = self.build_path()
        build_path.mkdir(parents=True, exist_ok=True)

        click.echo("Generating cmake files.")
        run_shell(f"cmake {def_str} {self.src}", cwd=build_path, extra_env=self.env)

    def build(self, targets: list[str], regenerate: bool = False) -> None:
        path = self.build_path()
        if regenerate:
            shutil.rmtree(path)
        if not path.exists():
            self._generate_cmake()
        run_shell(f"make -j$(($(nproc)-1)) {' '.join(targets)}", cwd=path, extra_env=self.env)
