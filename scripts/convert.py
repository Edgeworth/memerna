#!/usr/bin/env python3
import argparse

from common import read_file
from rna import RNA


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--to', choices=['db', 'ct'], required=True)
  parser.add_argument('filename')
  args = parser.parse_args()

  rna = RNA.from_any_file(read_file(args.filename))

  if args.to == 'ct':
    print(rna.to_ct_file())
  elif args.to == 'db':
    print(rna.to_db_file())


if __name__ == "__main__":
  main()
