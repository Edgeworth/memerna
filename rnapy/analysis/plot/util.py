from pathlib import Path
from typing import Any

from matplotlib import pyplot as plt
import seaborn as sns


def set_style() -> None:
    sns.set_theme(context="notebook", palette="deep")


def get_subplot_grid(
    n: int,
    sharex: bool = False,
    sharey: bool = False,
    inches: float = 3.0,
) -> tuple[plt.Figure, list[plt.Axes]]:
    SPLITTINGS = [(0, 0), (1, 1), (1, 2), (2, 2), (2, 2), (2, 3), (2, 3), (3, 3), (3, 3), (3, 3)]

    factors = SPLITTINGS[n]
    if factors == (1, 1):
        f, axes = plt.subplots(n, sharey=sharey, sharex=sharex)
    else:
        f, axes = plt.subplots(factors[0], factors[1], sharey=sharey, sharex=sharex)
        axes = axes.flatten()
    if n == 1:
        axes = [axes]
    f.tight_layout()

    # Hide subplots that are not used
    for i in range(n, factors[0] * factors[1]):
        axes[i].clear()
        axes[i].set_axis_off()
        axes[i].get_xaxis().set_visible(False)
        axes[i].get_yaxis().set_visible(False)

    f.set_size_inches(factors[1] * inches, factors[0] * inches)
    return f, axes


def get_marker(idx: int) -> dict[str, Any]:
    MARKERS = " ov^sp*+xD|"
    return {"marker": MARKERS[idx % len(MARKERS)], "markersize": 5, "markevery": 5}


def save_figure(f: plt.Figure, path: Path) -> None:
    f.tight_layout()
    f.savefig(path, dpi=300)
    plt.close(f)


def set_up_axis(ax: plt.Axes, names: tuple[str, str], legend: bool = True) -> None:
    for i, axis in [(0, ax.xaxis), (1, ax.yaxis)]:
        if names[i] is not None:
            axis.set_label_text(names[i])
    if legend:
        ax.legend(loc=0, fontsize="medium")


def set_up_figure(f: plt.Figure, names: tuple[str, str], legend: bool = True) -> None:
    f.suptitle(f"{names[0]} vs {names[1]}", y=1.00)
    for ax in f.get_axes():
        set_up_axis(ax, names, legend)
