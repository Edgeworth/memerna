# Copyright 2016 Eliot Courtney.
import os
import resource
import subprocess
import sys



def human_size(b, binary=True):
    def fmt(f):
        return (f"{f:.2f}").rstrip("0").rstrip(".")

    units = ["B", "KiB", "MiB", "GiB"]
    base = 1024
    if not binary:
        units = ["B", "KB", "MB", "GB"]
        base = 1000
    for unit in units[:-1]:
        if abs(b) < base:
            return f"{fmt(b)} {unit}"
        b /= base
    return f"{fmt(b)} {units[-1]}"
