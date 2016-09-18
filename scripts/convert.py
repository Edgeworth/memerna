#!/usr/bin/env python3
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
