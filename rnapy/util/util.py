# Copyright 2022 Eliot Courtney.
import hashlib
import inspect
import json
from typing import Any


def fn_args() -> dict[str, Any]:
    """Gets the arguments of the function that called this function."""
    frame = inspect.currentframe()
    if frame is None:
        raise RuntimeError("Cannot get frame.")
    if frame.f_back is None:
        raise RuntimeError("Cannot get prev.")
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


def stable_hash(val: Any) -> int:
    val = json.dumps(val, ensure_ascii=False, sort_keys=True, indent=None, separators=(",", ":"))
    val = hashlib.md5(val.encode("utf-8")).digest()
    return int.from_bytes(val, "big")
