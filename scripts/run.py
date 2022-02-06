#!/usr/bin/env python3
import cloup
from scripts.programs.analysis.compare_partition import compare_partition
from scripts.programs.convert_format import convert_format
from scripts.programs.crop_image import crop_image
from scripts.programs.harness import harness
from scripts.programs.parse_rnastructure_datatables import parse_rnastructure_datatables


CONTEXT_SETTINGS = cloup.Context.settings(
    show_constraints=True,
    show_subcommand_aliases=True,
    formatter_settings=cloup.HelpFormatter.settings(
        theme=cloup.HelpTheme.dark(),
    ),
)


@cloup.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


cli.section("Analysis", compare_partition, harness)
cli.section("Conversion", convert_format, crop_image, parse_rnastructure_datatables)

if __name__ == "__main__":
    cli()
