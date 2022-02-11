# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
from enum import Enum
from pathlib import Path


class Compiler(str, Enum):
    DEFAULT = "default"
    GCC = "gcc"
    CLANG = "clang"
    AFL_LTO = "afl-lto"
    AFL_FAST = "afl-fast"


class Sanitizer(str, Enum):
    NONE = "none"
    ASAN = "asan"
    TSAN = "tsan"
    UBSAN = "ubsan"


class BuildKind(str, Enum):
    DEBUG = "debug"
    RELEASE = "release"
    RELWITHDEBINFO = "relwithdebinfo"


@dataclass
class BuildCfg:
    src: Path
    prefix: Path
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
        return ident

    def build_path(self) -> Path:
        return self.prefix / "memerna" / self.ident()
