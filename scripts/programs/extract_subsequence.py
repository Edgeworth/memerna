#!/usr/bin/env python3
# Copyright 2016 Eliot Courtney.
import argparse

from scripts.common import read_file


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("starts", type=int, nargs="+")
    args = parser.parse_args()

    assert len(args.starts) % 2 == 0
    seq, pairs = read_file(args.filename).split("\n")
    for i in range(0, len(args.starts), 2):
        st, en = args.starts[i], args.starts[i + 1]
        print(f"Range [{int(st)}, {int(en)}]:\n  {seq[st:en + 1]}\n  {pairs[st:en + 1]}")
    print("Remember to reverse sequences sometimes.")


if __name__ == "__main__":
    main()
