# Copyright 2016 Eliot Courtney.
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from scripts.plot.load_data import colmap
from scripts.plot.plot_common import get_subplot_grid, set_up_figure, do_table


def subopt_distribution(all_ds, subopts):
  xid = 'strucpersec'
  frames = {}
  for subopt in subopts:
    ds_frames = [ds[subopt][ds[subopt][xid] != float('inf')] for ds in all_ds]
    frames[subopt] = pd.concat(ds_frames)

  do_table(frames, ['strucpersec'], True)
  f, axes = get_subplot_grid(len(frames))
  for i, frame_id in enumerate(sorted(frames.keys())):
    sns.distplot(frames[frame_id][xid], kde=False, bins=50, ax=axes[i], label=frame_id)

  f.tight_layout()
  for ax in axes:
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
  set_up_figure(f, names=(colmap[xid], None), legend=True)
  f.suptitle('Number of structures per second distribution', y=1.00)
  return f


def subopt_perf_results(ds):
  fmap_by_len = {name: frame.groupby('length') for name, frame in ds.fmap.items()}

  # savefig_local(
  #   ds.name, 'strucpersec',
  #   do_quantity_plot(fmap_by_len, 'length', 'strucpersec'))

  # savefig_local(
  #   ds.name, 'real',
  #   do_quantity_plot(fmap_by_len, 'length', 'real'))
  # savefig_local(
  #   ds.name, 'numstruc',
  #   do_quantity_plot(fmap_by_len, 'length', 'numstruc'))
  # savefig_local(
  #   ds.name, 'maxrss',
  #   do_quantity_plot(fmap_by_len, 'length', 'maxrss'))

  # savefig_local(
  #   ds.name, 'usersys',
  #   do_quantity_plot(fmap_by_len, 'length', ['real', 'usersys']))

  # f1, f2 = do_quantity_log_plot(fmap_by_len, 'length', 'real')
  # savefig_local(ds.name, 'real_loglog_scatter', f1)
  # savefig_local(ds.name, 'real_loglog_bestfit', f2)
  #
  # f1, f2 = do_quantity_log_plot(fmap_by_len, 'length', 'maxrss', logx=False)
  # savefig_local(ds.name, 'maxrss_logy_scatter', f1)
  # savefig_local(ds.name, 'maxrss_logy_bestfit', f2)
  #
  # f1, f2 = do_quantity_log_plot(fmap_by_len, 'length', 'maxrss')
  # savefig_local(ds.name, 'maxrss_loglog_scatter', f1)
  # savefig_local(ds.name, 'maxrss_loglog_bestfit', f2)
