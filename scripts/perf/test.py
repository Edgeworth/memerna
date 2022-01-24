#!/usr/bin/env python3
# Copyright 2016 E.
from scripts.common import *


def main():
    os.system("./scripts/run.py -b random -f -m")
    os.system("./scripts/run.py -b random_large -f -m")
    for delta in [1, 2, 3, 4, 5, 6, 10, 11, 12, 13]:
        os.system("./scripts/run.py -b random -s %d -m" % delta)
        os.system("./scripts/run.py -b random -s %d -ms" % delta)


if __name__ == "__main__":
    main()
