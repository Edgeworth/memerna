# Copyright 2022 Eliot Courtney.
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib as mpl

from rnapy.analysis.metrics import Dataset
from rnapy.analysis.plot.plots import Column, plot_mean_quantity
from rnapy.analysis.plot.util import save_figure
from rnapy.util.format import human_size


class FoldPerfPlotter:
    COLS: dict[str, Column] = {
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

    def _load_datasets(self) -> dict[str, Dataset]:
        datasets: dict[str, Dataset] = {}
        for path in self.input_dir.iterdir():
            dataset, program = path.stem.rsplit("_", maxsplit=1)
            df = pd.read_csv(path)
            datasets.setdefault(dataset, Dataset(dataset))
            datasets[dataset].dfs[program] = df
        return datasets

    def _path(self, ds: Dataset, name: str) -> Path:
        return self.output_dir / f"{ds.name}_{name}.png"

    def _plot(self, ds: Dataset) -> None:
        f = plot_mean_quantity(ds, self.COLS["length"], self.COLS["real_sec"])
        save_figure(f, self._path(ds, "test"))

    def run(self) -> None:
        datasets = self._load_datasets()
        for ds in datasets.values():
            self._plot(ds)
            break
