from dataclasses import dataclass

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from rnapy.analysis.metrics import Dataset
from rnapy.analysis.plot.util import get_marker
from rnapy.analysis.plot.util import get_subplot_grid
from rnapy.analysis.plot.util import set_up_figure
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf


@dataclass
class Column:
    idx: str
    name: str
    formatter: mpl.ticker.FuncFormatter | None = None

    def __str__(self) -> str:
        return self.idx


def plot_mean_quantity(ds: Dataset, xcol: Column, ycols: list[Column] | Column) -> plt.Figure:
    if not isinstance(ycols, list):
        ycols = [ycols]
    f, ax = plt.subplots(1)

    for idx, did in enumerate(sorted(ds.keys())):
        df = ds[did]

        for ycol in ycols:
            x, y = xcol.idx, ycol.idx
            df = df[[x, y]].groupby(x)
            with sns.color_palette(n_colors=len(ds)) as palette:
                sns.lineplot(df.mean(), x=x, y=y, label=did, ax=ax, **get_marker(idx))
                low, high = df[y].min(), df[y].max()
                ax.fill_between(list(sorted(df.groups)), low, high, alpha=0.2, color=palette.pop(0))

    set_up_figure(f, names=(xcol.name, ycols[0].name))
    if ycols[0].formatter:
        ax.yaxis.set_major_formatter(ycols[0].formatter)
    return f


def plot_mean_log_quantity(  # noqa: too-many-locals
    ds: Dataset,
    xcol: Column,
    ycol: Column,
    logx: bool = True,
    logy: bool = True,
) -> plt.Figure:
    EP = 1e-2

    f, axes = get_subplot_grid(len(ds), True, True)
    palette = sns.color_palette(n_colors=len(ds))
    for i, did in enumerate(sorted(ds.keys())):
        x, y = xcol.idx, ycol.idx

        df = ds[did][[x, y]].groupby(x).mean()
        df = df[df[y] > EP].reset_index()
        if logx:
            df[x] = df[x].apply(np.log10)
        if logy:
            df[y] = df[y].apply(np.log10)
        mod = smf.ols(f"{y} ~ {x}", data=df)
        res = mod.fit()

        b, a = res.params.tolist()
        sign = "-" if b < 0 else "+"
        label = f"{did}\n${a:.5f}x {sign} {abs(b):.2f}$\n$R^2 = {res.rsquared:.3f}$"
        sns.regplot(x=x, y=y, label=label, data=df, fit_reg=False, ax=axes[i], color=palette.pop(0))
        sm.graphics.abline_plot(model_results=res, ax=axes[i], c=(0, 0, 0, 0.8), **get_marker(i))

    names = [xcol.name, ycol.name]
    if logx:
        names[0] = f"log({names[0]})"
    if logy:
        names[1] = f"log({names[1]})"
    set_up_figure(f, names=(names[0], names[1]))

    return f
