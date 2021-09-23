#!/usr/bin/env python3
# Copyright 2021 Eliot Courtney.
import argparse
import decimal
import math

from scripts.common import read_file


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('file0', type=str)
  parser.add_argument('file1', type=str)
  args = parser.parse_args()

  p0 = [decimal.Decimal(i) for i in read_file(args.file0).split()]
  p1 = [decimal.Decimal(i) for i in read_file(args.file1).split()]
  rms = math.sqrt(sum([(p0[i] - p1[i])*(p0[i] - p1[i]) for i in range(len(p0))]) / len(p0))
  largest_diff = max(abs(p0[i] - p1[i]) for i in range(len(p0)))
  print('rms: %.20f' % rms)
  print('largest diff: %.20f' % largest_diff)


if __name__ == "__main__":
  main()
