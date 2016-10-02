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
import argparse

from common import read_file


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('filename')
  parser.add_argument('starts', type=int, nargs='+')
  args = parser.parse_args()

  assert (len(args.starts) % 2 == 0)
  seq, pairs = read_file(args.filename).split('\n')
  for i in range(0, len(args.starts), 2):
    st, en = args.starts[i], args.starts[i + 1]
    print('Range [%d, %d]:\n  %s\n  %s' % (st, en, seq[st:en + 1], pairs[st:en + 1]))
  print('Remember to reverse sequences sometimes.')


if __name__ == "__main__":
  main()
