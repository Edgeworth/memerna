# Copyright 2016 Eliot Courtney.
import pandas as pd
import scipy
from scripts.plot.load_data import colmap
from scripts.plot.plot_common import do_table
from scripts.plot.plot_common import get_subplot_grid
from scripts.plot.plot_common import latex_table
from scripts.plot.plot_common import savefig_local
from scripts.plot.plot_common import set_up_figure
import seaborn as sns

TEXT_LOC = (0.2, 0.9)


def do_accuracy_plot(frames, xid):
    f, axes = get_subplot_grid(len(frames))
    for i, frame_id in enumerate(sorted(frames.keys())):
        sns.distplot(frames[frame_id][xid], ax=axes[i], label=frame_id)
        axes[i].annotate(
            frame_id,
            xy=TEXT_LOC,
            xytext=TEXT_LOC,
            size=10,
            textcoords="axes fraction",
        )

    set_up_figure(f, names=(colmap[xid], None), legend=False)
    f.suptitle("F-Score distribution", y=1.00)
    return f


def do_comparison_plot(frames, aid, bid):
    f, axes = get_subplot_grid(len(frames))
    for i, frame_id in enumerate(sorted(frames.keys())):
        frame = frames[frame_id]
        sns.kdeplot(frame[aid], frame[bid], shade=True, ax=axes[i], label=frame_id)
        axes[i].annotate(
            frame_id,
            xy=TEXT_LOC,
            xytext=TEXT_LOC,
            size=10,
            textcoords="axes fraction",
        )
    set_up_figure(f, names=(colmap[aid], colmap[bid]))
    return f


def do_mfe_table(frames, xid):
    data = [i[xid] for i in frames.values()]
    data = pd.concat(data, axis=1)
    num_unique = data.apply(pd.Series.nunique, axis=1)
    num_unique = num_unique.value_counts().sort_index()
    num_unique /= num_unique.sum()
    table = [
        ["Number of MFEs"] + num_unique.index.tolist(),
        ["Proportion"] + [f"{i * 100:.2f}\\%" for i in num_unique.tolist()],
    ]
    print("MFE DISTRIBUTION TABLE:")
    print(latex_table(table))


def do_ttest(frames, xid):
    frame_keys = sorted(frames.keys())
    table = [["Package"] + frame_keys]
    for i, a_fid in enumerate(frame_keys):
        print(
            f"Package {a_fid} on {xid}: kurtosis: {frames[a_fid][xid].kurtosis():f}, skew: {frames[a_fid][xid].skew():f}",
        )
        table.append([a_fid])
        for j, b_fid in enumerate(frame_keys):
            if i == j:
                table[-1].append("---")
                continue
            frameA, frameB = frames[a_fid][xid], frames[b_fid][xid]
            if i < j:
                frameA, frameB = frameB, frameA
            t, p = scipy.stats.ttest_rel(frameA, frameB)
            entry = f"{t:.3f}"
            if i < j:
                entry = f"{p:.3f}"
            table[-1].append(entry)
    print("TTEST TABLE:")
    print(latex_table(table))


def fold_accuracy_results(ds):
    savefig_local(ds.name, "fscore", do_accuracy_plot(ds.fmap, "fscore"))
    savefig_local(ds.name, "comparison", do_comparison_plot(ds.fmap, "ppv", "sensitivity"))
    do_mfe_table(ds.fmap, "mfe")
    do_table(ds.fmap, ["fscore", "ppv", "sensitivity"])
    do_ttest(ds.fmap, "fscore")
