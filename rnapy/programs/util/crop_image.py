# Copyright 2016 E.
import os
from pathlib import Path

import cloup


@cloup.command(aliases=["crop"])
@cloup.argument(
    "files",
    type=cloup.Path(dir_okay=False, exists=True, writable=True, resolve_path=True, path_type=Path),
    nargs=-1,
)
def crop_image(files: list[Path]) -> None:
    for i in files:
        os.system(f"convert {i} -trim {i}")
