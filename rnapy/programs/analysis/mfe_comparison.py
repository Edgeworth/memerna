# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any
import cloup
from rnapy.bridge.args import bridge_options
from rnapy.data.args import memevault_options


@cloup.command()
@bridge_options
@memevault_options
def mfe_comparison(memevault_path: Path, **kwargs: Any) -> None:
    pass
