# Copyright 2022 Eliot Courtney.
import os


def fix_path(path):
    return os.path.abspath(os.path.expanduser(path))


def read_file(name):
    return open(name, encoding="utf-8").read()


def write_file(name, data):
    with open(name, "w") as f:
        f.write(data)
