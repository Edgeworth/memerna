#!/usr/bin/env python3
# Copyright 2021, E.
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
import argparse
import decimal
import math

from common import read_file


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('file0', type=str)
  parser.add_argument('file1', type=str)
  args = parser.parse_args()

  p0 = [decimal.Decimal(i) for i in read_file(args.file0).split()]
  p1 = [decimal.Decimal(i) for i in read_file(args.file0).split()]
  rms = math.sqrt(sum([abs(p0[i] - p1[i]) for i in range(len(p0))]))
  print(rms)


if __name__ == "__main__":
  main()
