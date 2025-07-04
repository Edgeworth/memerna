# Copyright 2022 E.
import copy
import enum
import hashlib
import inspect
import json
import tempfile
from enum import StrEnum
from pathlib import Path
from typing import IO, Any, cast

import click
import cloup
import polars as pl

from rnapy.util.command import run_cmd


def strict_merge(*dicts: dict) -> dict:
    merged = {}
    for curdict in dicts:
        for key, value in curdict.items():
            if key in merged:
                raise ValueError(f"Key '{key}' exists in both dictionaries.")
            merged[copy.deepcopy(key)] = copy.deepcopy(value)
    return merged


def row_by_key(json_path: Path, data_keys: dict) -> dict[str, Any] | None:
    """Checks if a row with the given keys exists in the JSON file."""
    if not json_path.exists():
        return None

    df = pl.read_ndjson(json_path)
    for key, value in data_keys.items():
        if key not in df.columns:
            return None
        df = df.filter(pl.col(key) == value)
    if df.is_empty():
        return None
    return cast(dict[str, Any], df.row(0, named=True))


def append_ndjson(path: Path, df: pl.DataFrame) -> None:
    """Appends a DataFrame as a new line in an NDJSON file."""
    with path.open(mode="a") as f:
        df.write_ndjson(f)


def append_csv(path: Path, df: pl.DataFrame) -> None:
    include_header = not path.exists() or path.stat().st_size == 0
    with path.open(mode="a") as f:
        df.write_csv(f, include_header=include_header)


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


def fast_linecount(path: Path) -> int:
    res = run_cmd("wc", "-l", str(path))
    count = int(res.stdout.strip().split()[0])
    return count


def resolve_path(path: Path | str) -> Path:
    return Path(path).expanduser().resolve()


def named_tmpfile(mode: str) -> IO[Any]:
    # If this runs out of space, use TMPDIR to change it.
    return tempfile.NamedTemporaryFile(mode)


class EnumChoice(cloup.Choice):
    """A custom Choice class for StrEnum types."""

    def normalize_choice(self, choice: StrEnum, ctx: click.Context | None) -> str:
        normed_value = choice.value if isinstance(choice, enum.Enum) else str(choice)

        if ctx is not None and ctx.token_normalize_func is not None:
            normed_value = ctx.token_normalize_func(normed_value)

        if not self.case_sensitive:
            normed_value = normed_value.casefold()

        return normed_value


def enum_choice(enum: type[StrEnum]) -> cloup.Choice:
    """Returns a list of choices for a StrEnum."""
    return EnumChoice(list(enum))
