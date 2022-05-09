#!/usr/bin/env python3
import logging

import click_log
import cloup
from rnapy.programs.afl_fuzz import afl_fuzz
from rnapy.programs.afl_fuzz_min import afl_fuzz_min
from rnapy.programs.analysis.compare_partition import compare_partition
from rnapy.programs.build import build
from rnapy.programs.convert_format import convert_format
from rnapy.programs.crop_image import crop_image
from rnapy.programs.design import design
from rnapy.programs.harness import harness
from rnapy.programs.parse_rnastructure_datatables import parse_rnastructure_datatables

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