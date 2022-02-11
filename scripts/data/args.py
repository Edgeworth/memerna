from pathlib import Path

import click
import cloup
from cloup.constraints import RequireAtLeast
from scripts.data.memevault import MemeVault
from scripts.model.parse import db_to_secondary
from scripts.model.parse import seq_to_primary
from scripts.model.rna import Rna


def validate_rna_file(ctx, param, value):
    if not value:
        return value
    return Rna.from_any_file(value.read())


def validate_seq(ctx, param, value):
    if not value:
        return value
    return seq_to_primary(value)


def validate_db(ctx, param, value):
    if not value:
        return value
    return db_to_secondary(value)


data_options = cloup.option_group(
    "Data options",
    cloup.option(
        "--memevault-path",
        envvar="MEMEVAULT",
        show_envvar=True,
        type=cloup.Path(dir_okay=False, exists=True, path_type=Path),
    ),
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
    **kwargs,
) -> Rna:
    if memevault_rna and not memevault_path:
        raise click.UsageError("--memevault-path is required when --memevault-rna is specified")
    if sum(1 for i in [memevault_rna, file_rna, primary or secondary] if i) != 1:
        raise click.UsageError(
            "Specify only one of --memevault-rna, --file-rna, --primary/--secondary",
        )
    if file_rna:
        return file_rna
    if memevault_rna:
        return MemeVault(memevault_path, "archiveii")[memevault_rna]
    if primary or secondary:
        return Rna("cmd", primary, secondary)
    raise click.UsageError("Specify one of --memevault-rna, --file-rna, --primary/--secondary")
