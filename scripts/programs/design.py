# Copyright 2022 E.
from pathlib import Path

import cloup

from scripts.design.run import run_design


@cloup.command()
@cloup.option(
    "--output-path",
    default=Path.home() / "bin" / "design",
    type=cloup.Path(
        dir_okay=True,
        file_okay=False,
        exists=False,
        writable=True,
        resolve_path=True,
        path_type=Path,
    ),
)
def design(output_path: Path) -> None:
    run_design(output_path)
