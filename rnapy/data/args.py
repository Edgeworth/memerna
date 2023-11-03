# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any, TextIO

import click
import cloup
from cloup.constraints import RequireAtLeast

from rnapy.data.memevault import MemeVault
from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.parse.sequence import db_to_secondary, seq_to_primary
from rnapy.model.rna import Rna


def validate_rna_file(
    _ctx: click.Context, _param: click.Parameter, value: TextIO | None
) -> Rna | None:
    if value is None:
        return None
    return RnaParser.from_any_file(value.read())


def validate_seq(_ctx: click.Context, _param: click.Parameter, value: str | None) -> str | None:
    if not value:
        return value
    return seq_to_primary(value)


def validate_db(
    _ctx: click.Context, _param: click.Parameter, value: str | None
) -> list[int] | None:
    if value is None:
        return None
    return db_to_secondary(value)


memevault_options = cloup.option(
    "--memevault-path",
    envvar="MEMEVAULT",
    show_envvar=True,
    type=cloup.Path(dir_okay=False, exists=True, path_type=Path),
)

data_options = cloup.option_group(
    "Data options",
    memevault_options,
    RequireAtLeast(1)(
        cloup.option("-mr", "--memevault-rna"),
        cloup.option("-fr", "--file-rna", type=cloup.File("r"), callback=validate_rna_file),
        cloup.option("-pr", "--primary", callback=validate_seq),
        cloup.option("-sr", "--secondary", callback=validate_db),
    ),
)


def rna_from_args(
    memevault_rna: str | None,
    memevault_path: Path | None,
    file_rna: Rna | None,
    primary: str | None,
    secondary: list[int] | None,
    **_kwargs: Any,
) -> Rna:
    if memevault_rna and not memevault_path:
        raise click.UsageError("--memevault-path is required when --memevault-rna is specified")
    if sum(1 for i in [memevault_rna, file_rna, primary or secondary] if i) != 1:
        raise click.UsageError(
            "Specify only one of --memevault-rna, --file-rna, --primary/--secondary"
        )
    if file_rna:
        return file_rna
    if memevault_rna:
        assert memevault_path is not None
        return MemeVault(memevault_path, "archiveii")[memevault_rna]
    if primary or secondary:
        return Rna("cmd", primary, secondary)
    raise click.UsageError("Specify one of --memevault-rna, --file-rna, --primary/--secondary")
