# Copyright 2022 Eliot Courtney.
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib as mpl

from rnapy.analysis.metrics import Dataset
from rnapy.analysis.plot.plots import Column
from rnapy.util.format import human_size


class FoldPerfPlotter:
    COLMAP: dict[str, Column] = {
        "name": Column(idx="name", name="Name"),
        "length": Column(idx="length", name="Length (nuc)"),
        "real_sec": Column(idx="real_sec", name="Wall time (s)"),
        "user_sec": Column(idx="user_sec", name="User time (s)"),
        "sys_sec": Column(idx="sys_sec", name="Sys time (s)"),
        "maxrss_bytes": Column(
            idx="maxrss_bytes",
            name="Maximum RSS (B)",
            formatter=mpl.ticker.FuncFormatter(lambda x, pos: human_size(x, False)),
        ),
    }
    input_dir: Path
    output_dir: Path

    def __init__(self, input_dir: Path, output_dir: Path) -> None:
        self.input_dir = input_dir
        self.output_dir = output_dir
        sns.set(color_codes=True)

    def load_datasets(self) -> dict[str, Dataset]:
        datasets: dict[str, Dataset] = {}
        for path in self.input_dir.iterdir():
            dataset, program = path.stem.split("_")
            df = pd.read_csv(path)
            datasets.setdefault(dataset, Dataset(dataset))
            datasets[dataset].dfs[program] = df
        return datasets

    def run(self) -> None:
        ds = self.load_datasets()
