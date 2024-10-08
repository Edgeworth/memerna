# Copyright 2022 E.
from pathlib import Path
from typing import ClassVar

import matplotlib as mpl
import pandas as pd

from rnapy.analysis.metrics import Dataset
from rnapy.analysis.plot.plots import Column, plot_mean_log_quantity, plot_mean_quantity
from rnapy.analysis.plot.util import save_figure, set_style
from rnapy.util.format import human_size


class SuboptPerfPlotter:
    COLS: ClassVar[dict[str, Column]] = {
        "name": Column(idx="name", name="Name"),
        "length": Column(idx="length", name="Length (nuc)"),
        "real_sec": Column(idx="real_sec", name="Wall time (s)"),
        "user_sec": Column(idx="user_sec", name="User time (s)"),
        "sys_sec": Column(idx="sys_sec", name="Sys time (s)"),
        "maxrss_bytes": Column(
            idx="maxrss_bytes",
            name="Maximum RSS (B)",
            formatter=mpl.ticker.FuncFormatter(lambda x, _: human_size(x, False)),
        ),
    }
    input_dir: Path
    output_dir: Path

    def __init__(self, input_dir: Path, output_dir: Path) -> None:
        self.input_dir = input_dir
        self.output_dir = output_dir
        set_style()

    def _load_datasets(self) -> dict[str, Dataset]:
        datasets: dict[str, Dataset] = {}
        for path in self.input_dir.glob("*.results"):
            dataset, program = path.stem.rsplit("_", maxsplit=1)
            df = pd.read_csv(path)
            datasets.setdefault(dataset, Dataset(dataset))
            datasets[dataset].dfs[program] = df
        return datasets

    def _path(self, ds: Dataset, name: str) -> Path:
        return self.output_dir / f"{ds.name}_{name}.png"

    def _plot_quantity(self, ds: Dataset, name: str = "") -> None:
        for y in ["real_sec", "maxrss_bytes"]:
            f = plot_mean_quantity(ds, self.COLS["length"], self.COLS[y])
            save_figure(f, self._path(ds, name + y))

    def run(self) -> None:
        datasets = self._load_datasets()
        # Plot quantities
        for ds in datasets.values():
            self._plot_quantity(ds)

            for y in ["real_sec", "maxrss_bytes"]:
                f = plot_mean_log_quantity(ds, self.COLS["length"], self.COLS[y])
                save_figure(f, self._path(ds, f"{y}_log"))
