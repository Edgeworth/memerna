# Copyright 2016 Eliot Courtney.
import pandas as pd

from scripts.common import fix_path
from scripts.common import read_file

colmap = {
    "name": "Name",
    "length": "Length (nuc)",
    "real": "Wall time (s)",
    "usersys": "User+sys time (s)",
    "maxrss": "Maximum RSS (B)",
    "fscore": "F-Score",
    "ppv": "PPV",
    "sensitivity": "Sensitivity",
    "mfe": "MFE",
    "numstruc": "Number of structures",
    "strucpersec": "Structures per second",
}


def read_general_frame(frame, subset, cols, samecols, avgcols):
    all_rows = []
    for name, group in frame.groupby("name"):
        assert len(group) == 5
        if name not in subset:
            continue
        # Assert that some columns are all the same
        for colname in samecols:
            if not (group[colname] == group[colname].iloc[0]).all():
                assert (group[colname] == group[colname].iloc[0]).all()
        rows = [[v for i, v in enumerate(group.iloc[0]) if i != 1] for _ in range(3)]
        # Remove outliers and average
        for colname in avgcols:
            col = list(group[colname])
            col.sort()
            col = col[1:-1]
            for i in range(3):
                rows[i][cols.index(colname) - 1] = col[i]
        all_rows.extend(rows)
    cleaned_frame = pd.DataFrame(all_rows, columns=cols[:1] + cols[2:])
    assert len(all_rows) // 3 == len(subset)
    return cleaned_frame


def read_fold_frame(filename, subset):
    cols = [
        "name",
        "run",
        "length",
        "real",
        "usersys",
        "maxrss",
        "fscore",
        "ppv",
        "sensitivity",
        "mfe",
    ]
    samecols = ["length", "fscore", "ppv", "sensitivity", "mfe"]
    avgcols = ["real", "usersys", "maxrss"]
    return read_general_frame(
        pd.read_csv(filename, delimiter=" ", header=None, names=cols),
        subset,
        cols,
        samecols,
        avgcols,
    )


def read_subopt_frame(filename, subset):
    cols = ["name", "run", "length", "real", "usersys", "maxrss", "numstruc", "strucpersec"]
    samecols = ["length"]
    avgcols = ["real", "usersys", "maxrss", "numstruc", "strucpersec"]
    frame = pd.read_csv(filename, delimiter=" ", header=None, names=cols)
    frame["strucpersec"] = frame["numstruc"] / frame["real"]
    # Have to put in NaNs for failed runs
    extra_rows = []
    for s in subset:
        if not frame[frame["name"] == s].empty:
            continue
        extra_rows += [
            [s, i, int(s[4:]), float("inf"), float("inf"), float("inf"), float("inf"), float("inf")]
            for i in range(5)
        ]
    extra_frame = pd.DataFrame(extra_rows, columns=cols)
    return read_general_frame(frame.append(extra_frame), subset, cols, samecols, avgcols)


class DataSet:
    def __init__(self, name, fmap):
        self.name = name
        self.fmap = fmap


def load_subset_file(subset_filename):
    return {i.strip() for i in read_file(fix_path(subset_filename)).strip().splitlines()}


def read_fold_dataset(name, filename_map, subset_filename):
    subset = load_subset_file(subset_filename)
    frames = {name: read_fold_frame(filename, subset) for name, filename in filename_map.items()}
    return DataSet(name, frames)


def read_subopt_dataset(name, filename_map, subset_filename):
    subset = load_subset_file(subset_filename)
    frames = {name: read_subopt_frame(filename, subset) for name, filename in filename_map.items()}
    return DataSet(name, frames)
