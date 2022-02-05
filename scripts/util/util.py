# Copyright 2022 Eliot Courtney.
import inspect


def fn_args():
    """Gets the arguments of the function that called this function."""
    keys, _, _, local = inspect.getargvalues(inspect.currentframe().f_back)
    args = {k: local[k] for k in keys}
    args.update(local.get("kwargs", {}))
    return args
