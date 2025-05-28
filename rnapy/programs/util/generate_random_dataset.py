# Copyright 2025 Eliot Courtney.
from pathlib import Path
from typing import Any

import click
import cloup

from rnapy.data.args import memevault_options
from rnapy.data.memevault import MemeVault


@cloup.command()
@memevault_options
@cloup.option(
    "--dataset-name",
    type=str,
    required=True,
    help="Name for the new dataset (table name in the database).",
)
@cloup.option("--use-range/--no-use-range", default=False, help="Use a range of lengths.")
@cloup.option("--start-length", type=int, help="Starting length for RNA sequences.")
@cloup.option("--end-length", type=int, help="Ending length for RNA sequences (inclusive).")
@cloup.option("--step-size", type=int, help="Step size to increment RNA length.")
@cloup.option(
    "--count-per-size",
    type=int,
    default=1,
    show_default=True,
    help="Number of random RNAs to generate for each size.",
)
@cloup.argument(
    "extra-lengths",
    type=int,
    nargs=-1,
    help="Additional lengths to generate, specified as integers.",
)
def generate_random_dataset(
    memevault_path: Path,
    dataset_name: str,
    use_range: bool,
    start_length: int | None,
    end_length: int | None,
    step_size: int | None,
    count_per_size: int,
    extra_lengths: list[int],
    **_kwargs: Any,
) -> None:
    memevault = MemeVault(memevault_path, dataset_name)

    lengths = set(extra_lengths)
    if use_range:
        if start_length is None or end_length is None or step_size is None:
            raise click.BadParameter("When using --use-range, specify start, end, and step size.")
        if start_length > end_length:
            raise click.BadParameter("start-length must be less than or equal to end-length.")
        click.echo(
            f"Generating {dataset_name} from {start_length} to {end_length}, "
            f"step {step_size}, {count_per_size} per size."
        )
        lengths.update(range(start_length, end_length + 1, step_size))

    for length in sorted(lengths):
        for i in range(count_per_size):
            name = f"random_len{length}_num{i + 1}"
            click.echo(f"Adding {name}")
            memevault.add_random(length, name=name)
