# Copyright 2023 Eliot Courtney.
from functools import reduce
from pathlib import Path

import matplotlib as mpl
import pandas as pd
from scipy.stats import ttest_rel
from rnapy.analysis.metrics import Dataset
from rnapy.analysis.plot.plots import Column
from rnapy.analysis.plot.plots import plot_mean_log_quantity
from rnapy.analysis.plot.plots import plot_mean_quantity
from rnapy.analysis.plot.util import save_figure
from rnapy.analysis.plot.util import set_style
from rnapy.util.format import human_size


class FoldAccuracyPlotter:
    COLS: dict[str, Column] = {
        "name": Column(idx="name", name="Name"),
        "family": Column(idx="family", name="Family"),
        "sensitivity": Column(idx="sensitivity", name="Sensitivity"),
        "ppv": Column(idx="ppv", name="Positive predictive value"),
        "f1": Column(idx="f1", name="F1 score"),
    }
    input_dir: Path
    output_dir: Path

    def __init__(self, input_dir: Path, output_dir: Path) -> None:
        self.input_dir = input_dir
        self.output_dir = output_dir
        set_style()

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

    def _plot_quantity(self, ds: Dataset) -> None:
        for y in ["real_sec", "maxrss_bytes"]:
            f = plot_mean_quantity(ds, self.COLS["length"], self.COLS[y])
            save_figure(f, self._path(ds, y))

    def _get_parent_rnas(self, df: pd.DataFrame) -> list[str]:
        domained = df[df["name"].str.contains("domain", case=False)]["name"].tolist()
        parents = set()
        for name in domained:
            parent = "_".join(name.split("_")[:-1])
            assert parent in df["name"].values
            parents.add(parent)

        return list(parents)

    def run(self) -> None:
        datasets = self._load_datasets()
        for ds in datasets.values():
            print(f"Dataset: {ds.name}")
            for did in ds:
                df = ds[did]
                parents = self._get_parent_rnas(df)
                df = df[~df["name"].isin(parents)]
                # # print summary of ds
                # # group by family
                gp = df.groupby("family")
                print(f"dataset {ds.name} program {did}")
                print(f"{gp['f1'].mean()}")

            print("paired t-tests:")
            df1 = ds["memerna-t04p2"].groupby("family")
            df2 = ds["memerna-t22p2"].groupby("family")
            for family, _ in df1:
                f1 = df1.get_group(family)["f1"]
                f2 = df2.get_group(family)["f1"]
                t_statistic, p_value = ttest_rel(f1, f2)
                print(f"Family {family}: t-statistic={t_statistic}, p-value={p_value}")
        # Plot quantities
        # for ds in datasets.values():
        #     self._plot_quantity(ds)
