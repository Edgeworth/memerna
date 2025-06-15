#!/usr/bin/env python3
import logging

import click_log
import cloup
from dotenv import load_dotenv

from rnapy.programs.afl_fuzz import afl_fuzz
from rnapy.programs.afl_fuzz_min import afl_fuzz_min
from rnapy.programs.analysis.fold_accuracy import run_fold_accuracy
from rnapy.programs.analysis.fold_perf import run_fold_perf
from rnapy.programs.analysis.subopt_perf import run_subopt_perf
from rnapy.programs.build import build
from rnapy.programs.harness import harness
from rnapy.programs.util.convert_format import convert_format
from rnapy.programs.util.generate_random_dataset import generate_random_dataset

CONTEXT_SETTINGS = cloup.Context.settings(
    show_constraints=True,
    show_subcommand_aliases=True,
    show_default=True,
    formatter_settings=cloup.HelpFormatter.settings(theme=cloup.HelpTheme.dark()),
)

load_dotenv()
logger = logging.getLogger()
click_log.basic_config(logger)


@cloup.group(context_settings=CONTEXT_SETTINGS)
@click_log.simple_verbosity_option(logger)
def cli() -> None:
    pass


cli.section("Runner", harness, run_fold_perf, run_fold_accuracy, run_subopt_perf)
cli.section("Build", build, afl_fuzz, afl_fuzz_min)
cli.section("Utility", convert_format, generate_random_dataset)

if __name__ == "__main__":
    cli()
