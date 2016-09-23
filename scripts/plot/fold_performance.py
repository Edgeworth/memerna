import matplotlib as mpl
import numpy as np
import seaborn as sns
import statsmodels as sm
import statsmodels.formula.api as smf
import statsmodels.graphics.regressionplots
from matplotlib import pyplot as plt

from common import human_size
from plot.load_data import colmap
from plot.plot_common import EP, set_up_figure, savefig_local, get_subplot_grid


def do_quantity_log_plot(frames, xid, yid, logx=True, logy=True):
  f, axes = get_subplot_grid(len(frames), True, True)
  f2, ax = plt.subplots(1)

  palette = sns.color_palette(n_colors=len(frames))
  for i, frame_id in enumerate(sorted(frames.keys())):
    data = frames[frame_id][[xid, yid]].mean()
    data = data[data[yid] > EP]
    if logx:
      data[xid] = data[xid].apply(np.log10)
    if logy:
      data[yid] = data[yid].apply(np.log10)
    mod = smf.ols('%s ~ %s' % (yid, xid), data=data)
    res = mod.fit()

    label = '%s\n$R^2 = %.3f$' % (frame_id, res.rsquared)
    sns.regplot(xid, yid, label=label, data=data, fit_reg=False, ax=axes[i])
    sm.graphics.regressionplots.abline_plot(
      model_results=res, ax=axes[i], c=(0, 0, 0, 0.8))

    eq_label = '{0}: ${2:.5f}x + {1:.2f}$'.format(frame_id, *res.params.tolist())
    sm.graphics.regressionplots.abline_plot(
      model_results=res, ax=ax, label=eq_label, c=palette[i])

  ax.set_xlim(axes[0].get_xlim())
  ax.set_ylim(axes[0].get_ylim())
  names = [colmap[xid], colmap[yid]]
  if logx:
    names[0] = 'log(%s)' % names[0]
  if logy:
    names[1] = 'log(%s)' % names[1]
  set_up_figure(f, names=names)
  set_up_figure(f2, names=names)

  return f, f2


def do_quantity_plot(frames, xid, yid):
  f, ax = plt.subplots(1)

  palette = sns.color_palette(n_colors=len(frames))
  for frame_id in sorted(frames.keys()):
    frame = frames[frame_id]
    frame[[xid, yid]].mean().plot(x=xid, y=yid, kind='line', ax=ax, label=frame_id)
    low, high = frame[yid].min(), frame[yid].max()
    ax.fill_between([i for i in sorted(frame.groups)], low, high, alpha=0.2, color=palette.pop(0))

  set_up_figure(f, names=(colmap[xid], colmap[yid]))
  if yid == 'maxrss':
    ax.yaxis.set_major_formatter(
      mpl.ticker.FuncFormatter(lambda x, pos: human_size(x, False)))
  return f


def fold_perf_results(ds):
  fmap_by_len = {name: frame.groupby('length') for name, frame in ds.fmap.items()}

  savefig_local(
    ds.name, 'real',
    do_quantity_plot(fmap_by_len, 'length', 'real'))
  savefig_local(
    ds.name, 'maxrss',
    do_quantity_plot(fmap_by_len, 'length', 'maxrss'))
  f1, f2 = do_quantity_log_plot(fmap_by_len, 'length', 'real')
  savefig_local(ds.name, 'real_loglog_scatter', f1)
  savefig_local(ds.name, 'real_loglog_bestfit', f2)

  f1, f2 = do_quantity_log_plot(fmap_by_len, 'length', 'maxrss', logx=False)
  savefig_local(ds.name, 'maxrss_logy_scatter', f1)
  savefig_local(ds.name, 'maxrss_logy_bestfit', f2)

  f1, f2 = do_quantity_log_plot(fmap_by_len, 'length', 'maxrss')
  savefig_local(ds.name, 'maxrss_loglog_scatter', f1)
  savefig_local(ds.name, 'maxrss_loglog_bestfit', f2)
