#!/usr/bin/env python3
# Copyright 2016 Eliot Courtney.
from scripts.common import *


def main():
  os.system('./scripts/run.py -b random -f -k')
  os.system('./scripts/run.py -b random_large -f -k')
  for delta in [1, 2, 3, 4, 5, 6, 10, 11, 12, 13]:
    os.system('./scripts/run.py -b random -s %d -k' % delta)
    os.system('./scripts/run.py -b random -s %d -ks' % delta)


if __name__ == '__main__':
  main()
