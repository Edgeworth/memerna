# Copyright 2021 Eliot Courtney.
import decimal
import math

import click
import cloup


@cloup.command()
@cloup.argument("p0", type=cloup.File("r"))
@cloup.argument("p1", type=cloup.File("r"))
def compare_partition(p0_str: str, p1_str: str) -> None:
    p0 = [decimal.Decimal(i) for i in p0_str.split()]
    p1 = [decimal.Decimal(i) for i in p1_str.split()]
    rms = math.sqrt(sum((p0[i] - p1[i]) * (p0[i] - p1[i]) for i in range(len(p0))) / len(p0))
    largest_diff = max(abs(p0[i] - p1[i]) for i in range(len(p0)))
    click.echo(f"rms: {rms:.20f}")
    click.echo(f"largest diff: {largest_diff:.20f}")
