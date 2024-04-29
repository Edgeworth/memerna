from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import matplotlib as mpl
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from bidict import bidict
from matplotlib import pyplot as plt

from rnapy.analysis.metrics import Dataset
from rnapy.analysis.plot.util import get_marker, get_subplot_grid, set_up_figure
from rnapy.util.util import stable_hash


@dataclass
class Column:
    idx: str
    name: str
    formatter: mpl.ticker.FuncFormatter | None = None

    def __str__(self) -> str:
        return self.idx


def _color(name: str, palette: Sequence[Any] | None = None) -> Any:
    color_map: bidict = getattr(_color, "color_map", bidict())
    _color.color_map = color_map  # type: ignore[attr-defined]

    print(color_map)

    if palette is None:
        palette = sns.utils.get_color_cycle()

    if name not in color_map:
        start_idx = stable_hash(name) % len(palette)
        idx = start_idx
        while palette[idx] in color_map.inverse:
            idx = (idx + 1) % len(palette)
            if idx == start_idx:
                raise ValueError("No free color available in the palette.")

        color_map[name] = palette[idx]

    return color_map[name]


def plot_mean_quantity(ds: Dataset, xcol: Column, ycols: list[Column] | Column) -> plt.Figure:
    if not isinstance(ycols, list):
        ycols = [ycols]
    f, ax = plt.subplots(1)

    for idx, did in enumerate(sorted(ds.keys())):
        df = ds[did]

        for ycol in ycols:
            x, y = xcol.idx, ycol.idx
            df = df[[x, y]].groupby(x)
            sns.lineplot(
                df.mean(), x=x, y=y, label=did, ax=ax, color=_color(did), **get_marker(idx)
            )
            low, high = df[y].min(), df[y].max()
            ax.fill_between(sorted(df.groups), low, high, alpha=0.2, color=_color(did))

    set_up_figure(f, names=(xcol.name, ycols[0].name))
    if ycols[0].formatter:
        ax.yaxis.set_major_formatter(ycols[0].formatter)
    ax.legend(loc="best", framealpha=0.5)
    return f


def plot_mean_log_quantity(
    ds: Dataset, xcol: Column, ycol: Column, logx: bool = True, logy: bool = True
) -> plt.Figure:
    ep = 1e-2

    f, axes = get_subplot_grid(len(ds), sharex=True, sharey=True)
    for i, did in enumerate(sorted(ds.keys())):
        x, y = xcol.idx, ycol.idx

        df = ds[did][[x, y]].groupby(x).mean()
        df = df[df[y] > ep].reset_index()
        if logx:
            df[x] = df[x].apply(np.log10)
        if logy:
            df[y] = df[y].apply(np.log10)
        mod = smf.ols(f"{y} ~ {x}", data=df)
        res = mod.fit()

        b, a = res.params.tolist()
        sign = "-" if b < 0 else "+"
        label = f"{did}\n${a:.5f}x {sign} {abs(b):.2f}$\n$R^2 = {res.rsquared:.3f}$"
        sns.regplot(x=x, y=y, label=label, data=df, fit_reg=False, ax=axes[i], color=_color(did))
        sm.graphics.abline_plot(
            model_results=res, ax=axes[i], c=(0, 0, 0, 0.8), color=_color(did), **get_marker(i)
        )

    names = [xcol.name, ycol.name]
    if logx:
        names[0] = f"log({names[0]})"
    if logy:
        names[1] = f"log({names[1]})"
    set_up_figure(f, names=(names[0], names[1]))

    return f
