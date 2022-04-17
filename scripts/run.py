#!/usr/bin/env python3
import logging

import click_log
import cloup
from scripts.programs.afl_fuzz import afl_fuzz
from scripts.programs.afl_fuzz_min import afl_fuzz_min
from scripts.programs.analysis.compare_partition import compare_partition
from scripts.programs.build import build
from scripts.programs.convert_format import convert_format
from scripts.programs.crop_image import crop_image
from scripts.programs.design import design
from scripts.programs.harness import harness
from scripts.programs.parse_rnastructure_datatables import parse_rnastructure_datatables

CONTEXT_SETTINGS = cloup.Context.settings(
    show_constraints=True,
    show_subcommand_aliases=True,
    formatter_settings=cloup.HelpFormatter.settings(
        theme=cloup.HelpTheme.dark(),
    ),
)

logger = logging.getLogger()
click_log.basic_config(logger)


@cloup.group(context_settings=CONTEXT_SETTINGS)
@click_log.simple_verbosity_option(logger)
def cli() -> None:
    pass


cli.section("Analysis", compare_partition, harness)
cli.section("Build", build, afl_fuzz, afl_fuzz_min)
cli.section("Conversion", convert_format, crop_image, parse_rnastructure_datatables)
cli.section("Design", design)

if __name__ == "__main__":
    cli()
