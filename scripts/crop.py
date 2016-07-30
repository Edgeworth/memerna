#!/usr/bin/env python3
import argparse
import os


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('files', nargs='+')
  args = parser.parse_args()
  for i in args.files:
    os.system('convert %s -trim %s' % (i, i))


if __name__ == '__main__':
  main()
