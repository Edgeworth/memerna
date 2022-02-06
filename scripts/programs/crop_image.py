# Copyright 2016 Eliot Courtney.
import os

import cloup


@cloup.command(aliases=["crop"])
@cloup.argument("files", type=cloup.File("r"), nargs=-1)
def crop_image(files):
    for i in files:
        os.system(f"convert {i} -trim {i}")
