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
