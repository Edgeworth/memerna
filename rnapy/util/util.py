# Copyright 2022 Eliot Courtney.
import inspect
from typing import Any


def fn_args() -> dict[str, Any]:
    """Gets the arguments of the function that called this function."""
    frame = inspect.currentframe()
    assert frame is not None and frame.f_back is not None
    keys, _, _, local = inspect.getargvalues(frame.f_back)
    args = {k: local[k] for k in keys}
    args.update(local.get("kwargs", {}))
    args.update(local.get("_kwargs", {}))
    return args


def human_size(num_bytes: int, binary: bool = True) -> str:
    def float_fmt(f: float) -> str:
        return f"{f:.2f}".rstrip("0").rstrip(".")

    units = ["B", "KiB", "MiB", "GiB"]
    base = 1024
    if not binary:
        units = ["B", "KB", "MB", "GB"]
        base = 1000

    v = float(num_bytes)
    for unit in units[:-1]:
        if abs(v) < base:
            return f"{float_fmt(v)} {unit}"
        v /= base
    return f"{float_fmt(v)} {units[-1]}"
