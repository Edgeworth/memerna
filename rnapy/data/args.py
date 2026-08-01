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
    ctx: click.Context, param: click.Parameter, value: TextIO | None
) -> Rna | None:
    if value is None:
        return None
    try:
        return RnaParser.from_any_file(value.read())
    except ValueError as exc:
        raise click.BadParameter(str(exc), ctx=ctx, param=param) from exc


def validate_seq(ctx: click.Context, param: click.Parameter, value: str | None) -> str | None:
    if not value:
        return value
    try:
        return seq_to_primary(value)
    except ValueError as exc:
        raise click.BadParameter(str(exc), ctx=ctx, param=param) from exc


def validate_db(ctx: click.Context, param: click.Parameter, value: str | None) -> list[int] | None:
    if value is None:
        return None
    try:
        return db_to_secondary(value)
    except ValueError as exc:
        raise click.BadParameter(str(exc), ctx=ctx, param=param) from exc


def memevault_path_option(*, required: bool, exists: bool = True) -> Any:
    return cloup.option(
        "--memevault-path",
        envvar="MEMEVAULT",
        show_envvar=True,
        required=required,
        type=cloup.Path(dir_okay=False, exists=exists, path_type=Path),
    )


memevault_options = cloup.option_group(
    "memevault options",
    memevault_path_option(required=False),
    cloup.option("--memevault-dataset", default="archiveii", show_default=True),
)


def data_options(fn: Any) -> Any:
    return memevault_options(
        cloup.option_group(
            "Data options",
            RequireAtLeast(1)(
                cloup.option("-mr", "--memevault-rna"),
                cloup.option("-fr", "--file-rna", type=cloup.File("r"), callback=validate_rna_file),
                cloup.option("-pr", "--primary", callback=validate_seq),
                cloup.option("-sr", "--secondary", callback=validate_db),
            ),
        )(fn)
    )


def rna_from_args(
    memevault_rna: str | None,
    memevault_path: Path | None,
    memevault_dataset: str,
    file_rna: Rna | None,
    primary: str | None,
    secondary: list[int] | None,
    **_kwargs: Any,
) -> Rna:
    if sum(1 for i in [memevault_rna, file_rna, primary or secondary] if i) != 1:
        raise click.UsageError(
            "Specify only one of --memevault-rna, --file-rna, --primary/--secondary"
        )
    if file_rna:
        return file_rna
    if memevault_rna:
        if memevault_path is None:
            raise click.UsageError("--memevault-path is required when --memevault-rna is specified")
        vault = MemeVault(memevault_path, memevault_dataset)
        if not vault.dataset_exists():
            raise click.BadParameter(
                f"dataset {memevault_dataset!r} does not exist", param_hint="--memevault-dataset"
            )
        try:
            return vault[memevault_rna]
        except KeyError as exc:
            raise click.BadParameter(
                f"RNA {memevault_rna!r} does not exist in dataset {memevault_dataset!r}",
                param_hint="--memevault-rna",
            ) from exc
    if primary or secondary:
        return Rna("cmd", primary, secondary)
    raise click.UsageError("Specify one of --memevault-rna, --file-rna, --primary/--secondary")
