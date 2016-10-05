#!/usr/bin/env python3
# Copyright 2016, E.
#
# This file is part of memerna.
#
# memerna is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with memerna.
# If not, see <http://www.gnu.org/licenses/>.
from common import *


def main():
  os.system('./scripts/run.py -b random -f -k')
  os.system('./scripts/run.py -b random_large -f -k')
  for delta in [1, 2, 3, 4, 5, 6, 10, 11, 12, 13]:
    # os.system('./scripts/run.py -b random -s %d -rh' % delta)
    # os.system('./scripts/run.py -b random -s %d -rd' % delta)
    # os.system('./scripts/run.py -b random -s %d -vd2' % delta)
    # os.system('./scripts/run.py -b random -s %d -vd3' % delta)
    # os.system('./scripts/run.py -b random -s %d -vd2s' % delta)
    # os.system('./scripts/run.py -b random -s %d -vd3s' % delta)
    # os.system('./scripts/run.py -b random -s %d -sjsmpi' % delta)
    # os.system('./scripts/run.py -b random -s %d -sjsmpis' % delta)
    os.system('./scripts/run.py -b random -s %d -k' % delta)
    # pass
  # TODO do next:

  # os.system('./scripts/run.py -b random -s %d -k' % delta)


if __name__ == '__main__':
  main()
