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
from rna import RNA


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--to', choices=['db', 'ct'], required=True)
  parser.add_argument('-f', '--filename', type=str)
  parser.add_argument('data', type=str, nargs='2', required=False)
  args = parser.parse_args()

  if args.filename:
    rna = RNA.from_any_file(read_file(args.filename))
  else:
    rna = RNA.from_name_seq_db('user', *args.data)
  if args.to == 'ct':
    print(rna.to_ct_file())
  elif args.to == 'db':
    print(rna.to_db_file())


if __name__ == "__main__":
  main()
