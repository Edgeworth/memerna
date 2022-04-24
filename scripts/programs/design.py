# Copyright 2022 Eliot Courtney.
from pathlib import Path

import cloup
from scripts.design.run import run


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
@cloup.option(
    "--checkpoint-path",
    required=False,
    type=cloup.Path(
        dir_okay=False,
        file_okay=True,
        exists=True,
        writable=False,
        resolve_path=True,
        path_type=Path,
    ),
)
def design(output_path: Path, checkpoint_path: Path | None) -> None:
    run(output_path, checkpoint_path)
