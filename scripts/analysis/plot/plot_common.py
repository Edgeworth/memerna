# Copyright 2016 Eliot Courtney.
import matplotlib as mpl
import numpy as np
import seaborn as sns
import statsmodels as sm
import statsmodels.formula.api as smf
import statsmodels.graphics.regressionplots
from matplotlib import pyplot as plt

from scripts.common import human_size
from scripts.plot.load_data import colmap

EP = 1e-2


def set_up_axis(ax, names, legend):
    for i, axis in [(0, ax.xaxis), (1, ax.yaxis)]:
        if names[i] is not None:
            axis.set_label_text(f"{names[i]}")
    if legend:
        ax.legend(loc=0, fontsize="medium")


def set_up_figure(f, names=None, legend=True):
    f.suptitle("%s vs %s" % tuple(reversed(names)), y=1.00)
    for ax in f.get_axes():
        set_up_axis(ax, names, legend)


def save_figure(f, name):
    f.tight_layout()
    f.savefig(name, dpi=150)


def savefig_local(dataset, name, f):
    save_figure(f, f"./build/benchmark_figures/{dataset}_{name}.png")
    plt.close(f)


def latex_table(rows):
    result = ""
    for row in rows:
        row = np.array(row).tolist()
        result += f"{' & '.join([str(i) for i in row])} \\\\\n"
    return result


SPLITTINGS = [(0, 0), (1, 1), (1, 2), (2, 2), (2, 2), (2, 3), (2, 3), (3, 3), (3, 3), (3, 3)]


def get_subplot_grid(n, sharex=False, sharey=False):
    factors = SPLITTINGS[n]
    if factors == (1, 1):
        f, axes = plt.subplots(n, sharey=sharey, sharex=sharex)
    else:
        f, axes = plt.subplots(factors[0], factors[1], sharey=sharey, sharex=sharex)
        axes = axes.flatten()
    if n == 1:
        axes = [axes]
    f.tight_layout()
    # if factors[0] * factors[1] != n:
    #   axes[-1].clear()
    #   axes[-1].set_axis_off()
    #   axes[-1].get_xaxis().set_visible(False)
    #   axes[-1].get_yaxis().set_visible(False)
    f.set_size_inches(factors[1] * 3, factors[0] * 3)
    return f, axes


def get_marker(idx):
    s = " ov^sp*+xD|"
    return {"marker": s[idx % len(s)], "markersize": 5, "markevery": 5}


def do_quantity_log_plot(frames, xid, yid, logx=True, logy=True):
    f, axes = get_subplot_grid(len(frames), True, True)

    palette = sns.color_palette(n_colors=len(frames))
    for i, frame_id in enumerate(sorted(frames.keys())):
        data = frames[frame_id][[xid, yid]].mean()
        data = data[data[yid] > EP]
        if logx:
            data[xid] = data[xid].apply(np.log10)
        if logy:
            data[yid] = data[yid].apply(np.log10)
        mod = smf.ols(f"{yid} ~ {xid}", data=data)
        res = mod.fit()

        label = "{0}\n${3:.5f}x + {2:.2f}$\n$R^2 = {1:.3f}$".format(
            frame_id, res.rsquared, *res.params.tolist()
        )
        sns.regplot(xid, yid, label=label, data=data, fit_reg=False, ax=axes[i])
        sm.graphics.regressionplots.abline_plot(
            model_results=res, ax=axes[i], c=(0, 0, 0, 0.8), **get_marker(i)
        )

    names = [colmap[xid], colmap[yid]]
    if logx:
        names[0] = f"log({names[0]})"
    if logy:
        names[1] = f"log({names[1]})"
    set_up_figure(f, names=names)

    return f


def do_quantity_plot(frames, xid, yids):
    if not isinstance(yids, list):
        yids = [yids]
    f, ax = plt.subplots(1)

    for idx, frame_id in enumerate(sorted(frames.keys())):
        frame = frames[frame_id]
        for yid in yids:
            palette = sns.color_palette(n_colors=len(frames))
            frame[[xid, yid]].mean().plot(
                x=xid, y=yid, kind="line", ax=ax, label=frame_id, **get_marker(idx)
            )
            low, high = frame[yid].min(), frame[yid].max()
            ax.fill_between(
                [i for i in sorted(frame.groups)], low, high, alpha=0.2, color=palette.pop(0)
            )

    set_up_figure(f, names=(colmap[xid], colmap[yids[0]]))
    if yids[0] == "maxrss":
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: human_size(x, False)))
    return f


def do_table(frames, ids, do_median=False):
    cols = ["mean", "SD"]
    if do_median:
        cols.append("median")
    table = [["Package"] + [f"{colmap[i]} {t}" for i in ids for t in cols]]
    for frame_id in sorted(frames.keys()):
        frame = frames[frame_id]
        for i in ids:
            cols = [frame[i].mean(), frame[i].std()]
            if do_median:
                cols.append(frame[i].median())
            table.append([frame_id] + [f"{val:.5f}" for val in cols])
    print("TABLE:")
    print(latex_table(table))
