#!/usr/bin/env python3
# Copyright 2016 Eliot Courtney.
import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()
    for i in args.files:
        os.system(f"convert {i} -trim {i}")


if __name__ == "__main__":
    main()
