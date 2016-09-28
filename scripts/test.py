#!/usr/bin/env python3
from common import *


def main():
  # DONE: 1, 2, 3, 6, 10
  for delta in [4, 5, 11, 12, 13]:
    os.system('./scripts/run.py -b random -s %d -rh' % delta)
    os.system('./scripts/run.py -b random -s %d -rd' % delta)
    os.system('./scripts/run.py -b random -s %d -vd2' % delta)
    os.system('./scripts/run.py -b random -s %d -vd3' % delta)
    os.system('./scripts/run.py -b random -s %d -vd2s' % delta)
    os.system('./scripts/run.py -b random -s %d -vd3s' % delta)

    #os.system('./scripts/run.py -b random -s %d -k' % delta)


if __name__ == '__main__':
  main()
