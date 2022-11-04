from dataclasses import dataclass
from matplotlib import pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns
import statsmodels as sm
import statsmodels.formula.api as smf
from rnapy.analysis.metrics import Dataset
from rnapy.analysis.plot.util import get_marker, get_subplot_grid, set_up_figure


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
                df.mean().plot(kind="line", ax=ax, label=did, **get_marker(idx))
                low, high = df[y].min(), df[y].max()
                ax.fill_between(list(sorted(df.groups)), low, high, alpha=0.2, color=palette.pop(0))

    set_up_figure(f, names=(xcol.name, ycols[0].name))
    if ycols[0].formatter:
        ax.yaxis.set_major_formatter(ycols[0].formatter)
    return f


def plot_quantity_log(
    ds: Dataset, xcol: Column, ycol: Column, logx: bool = True, logy: bool = True
) -> plt.Figure:
    EP = 1e-2

    f, axes = get_subplot_grid(len(ds), True, True)
    palette = sns.color_palette(n_colors=len(ds))
    for i, did in enumerate(sorted(ds.keys())):
        df = ds[did][[xcol, ycol]].mean()
        df = df[df[ycol] > EP]
        if logx:
            df[xcol] = df[xcol].apply(np.log10)
        if logy:
            df[ycol] = df[ycol].apply(np.log10)
        mod = smf.ols(f"{ycol} ~ {xcol}", data=df)
        res = mod.fit()

        label = "{0}\n${3:.5f}x + {2:.2f}$\n$R^2 = {1:.3f}$".format(
            did, res.rsquared, *res.params.tolist()
        )
        sns.regplot(xcol, ycol, label=label, data=df, fit_reg=False, ax=axes[i])
        sm.graphics.regressionplots.abline_plot(
            model_results=res, ax=axes[i], c=(0, 0, 0, 0.8), **get_marker(i)
        )

    names = [xcol.name, ycol.name]
    if logx:
        names[0] = f"log({names[0]})"
    if logy:
        names[1] = f"log({names[1]})"
    set_up_figure(f, names=(names[0], names[1]))

    return f
