# Copyright 2016 Eliot Courtney.
from scripts.build.config import BuildCfg
from scripts.build.config import Compiler
from scripts.build.config import Sanitizer
from scripts.util.command import run_shell


def get_compilers(compiler: Compiler) -> tuple[str, str]:
    if compiler == Compiler.DEFAULT:
        return "cc", "c++"
    if compiler == Compiler.GCC:
        return "gcc", "g++"
    if compiler == Compiler.CLANG:
        return "clang", "clang++"
    if compiler == Compiler.AFL_LTO:
        return "afl-clang-lto", "afl-clang-lto++"
    if compiler == Compiler.AFL_FAST:
        return "afl-clang-fast", "afl-clang-fast++"
    raise ValueError(f"Unknown compiler {compiler}")


def get_sanitizers(sanitizer: Sanitizer) -> dict[str, str]:
    if sanitizer == Sanitizer.NONE:
        return {}
    if sanitizer == Sanitizer.ASAN:
        return {"SANITIZE_ADDRESS": "ON"}
    if sanitizer == Sanitizer.TSAN:
        return {"SANITIZE_THREAD": "ON"}
    if sanitizer == Sanitizer.UBSAN:
        return {"SANITIZE_UNDEFINED": "ON"}
    raise ValueError(f"Unknown sanitizer {sanitizer}")


def generate_cmake(cfg: BuildCfg) -> None:
    cc, cxx = get_compilers(cfg.compiler)

    defs = {
        "CMAKE_C_COMPILER": cc,
        "CMAKE_CXX_COMPILER": cxx,
        "CMAKE_BUILD_TYPE": cfg.kind,
        "USE_RNASTRUCTURE": "ON" if cfg.rnastructure else "OFF",
        "USE_MPFR": "ON" if cfg.mpfr else "OFF",
        "USE_IWYU": "ON" if cfg.iwyu else "OFF",
        "FLOAT_BITS": f"{cfg.float_bits}",
    }
    defs.update(get_sanitizers(cfg.sanitizer))
    def_str = " ".join(f"-D {i}={k}" for i, k in defs.items())

    build_path = cfg.build_path()
    build_path.mkdir(parents=True, exist_ok=True)
    run_shell(f"cmake {def_str} {cfg.src}", cwd=build_path)
