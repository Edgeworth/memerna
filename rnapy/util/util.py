# Copyright 2022 Eliot Courtney.
import copy
import enum
import hashlib
import inspect
import json
import subprocess
import tempfile
from collections.abc import Callable
from concurrent.futures import Future, ThreadPoolExecutor, wait
from enum import StrEnum
from pathlib import Path
from typing import IO, Any

import click
import cloup
import polars as pl


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

    ndjson = pl.read_ndjson(json_path)
    for key, value in data_keys.items():
        if key not in ndjson.columns:
            return None
        if value is None:
            ndjson = ndjson.filter(pl.col(key).is_null())
        else:
            ndjson = ndjson.filter(pl.col(key) == value)
    if ndjson.is_empty():
        return None
    return ndjson.row(0, named=True)


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


def stable_hash(val: Any) -> int:
    val = json.dumps(val, ensure_ascii=False, sort_keys=True, indent=None, separators=(",", ":"))
    val = hashlib.md5(val.encode("utf-8")).digest()
    return int.from_bytes(val, "big")


def fast_linecount(path: Path, compressed: bool = False) -> int:
    cat_cmd = ["zstd", "-dc", str(path)] if compressed else ["cat", str(path)]
    cat_proc = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    wc_cmd = ["wc", "-l"]
    wc_proc = subprocess.Popen(wc_cmd, stdin=cat_proc.stdout, stdout=subprocess.PIPE)
    assert cat_proc.stdout is not None
    cat_proc.stdout.close()
    wc_output, _ = wc_proc.communicate()
    _, cat_stderr = cat_proc.communicate()
    if cat_proc.returncode != 0:
        raise RuntimeError(f"{cat_cmd[0]} failed: {cat_stderr.decode()}")
    return int(wc_output.strip())


def read_compressed(path: Path) -> str:
    result = subprocess.run(
        ["zstd", "-dc", str(path)],  # noqa: S607
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"zstd -dc failed: {result.stderr}")
    return result.stdout


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


def parallel_map(fn: Callable[..., Any], jobs: list[tuple[Any, ...]], max_workers: int) -> None:
    """Run fn(*args) for each args in jobs, using a thread pool if max_workers > 1.

    Handles Ctrl-C cleanly by cancelling pending futures."""
    if max_workers <= 1:
        for args in jobs:
            fn(*args)
        return

    pool = ThreadPoolExecutor(max_workers=max_workers)
    futures: list[Future[Any]] = [pool.submit(fn, *args) for args in jobs]
    try:
        remaining = set(futures)
        while remaining:
            done, remaining = wait(remaining, timeout=1)
            for f in done:
                f.result()
    except BaseException:
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    pool.shutdown()
