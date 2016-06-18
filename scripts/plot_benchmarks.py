#!/usr/bin/env python3
import matplotlib as mpl
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import statsmodels as sm
import statsmodels.formula.api as smf
import statsmodels.graphics.regressionplots
import statsmodels.stats.weightstats
import statsmodels.stats.diagnostic
from common import human_size
from matplotlib import pyplot as plt

sns.set(color_codes=True)

EP = 1e-2

cols = {
  'length': 'Length (nuc)',
  'real': 'Wall time (s)',
  'usersys': 'User+sys time (s)',
  'maxrss': 'Maximum RSS (B)',
  'fscore': 'F-Score',
  'ppv': 'PPV',
  'sensitivity': 'Sensitivity',
  'mfe': 'MFE'
}


def set_up_axis(ax, names, legend):
  for i, axis in [(0, ax.xaxis), (1, ax.yaxis)]:
    if names[i] is not None:
      axis.set_label_text('%s' % names[i])
  if legend:
    ax.legend(loc=0, fontsize='medium')


def set_up_figure(f, names=None, legend=True):
  f.suptitle('%s vs %s' % tuple(reversed(names)))
  for ax in f.get_axes():
    set_up_axis(ax, names, legend)


def save_figure(f, name):
  f.tight_layout()
  f.savefig(name)


def do_quantity_log_plot(frames, xid, yid, logx=True, logy=True):
  f, axes = plt.subplots(len(frames) // 2, len(frames) // 2, sharey=True, sharex=True)
  f2, ax = plt.subplots(1)
  axes = axes.flatten()

  for i, frame_id in enumerate(frames):
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
      model_results=res, ax=ax, label=eq_label, c=sns.color_palette()[i])

  ax.set_xlim(axes[0].get_xlim())
  ax.set_ylim(axes[0].get_ylim())
  names = [cols[xid], cols[yid]]
  if logx:
    names[0] = 'log(%s)' % names[0]
  if logy:
    names[1] = 'log(%s)' % names[1]
  set_up_figure(f, names=names)
  set_up_figure(f2, names=names)

  return f, f2


def do_quantity_plot(frames, xid, yid):
  f, ax = plt.subplots(1)

  for frame_id, frame in frames.items():
    frame[[xid, yid]].mean().plot(x=xid, y=yid, kind='line', ax=ax, label=frame_id)

  set_up_figure(f, names=(cols[xid], cols[yid]))
  if yid == 'maxrss':
    ax.yaxis.set_major_formatter(
      mpl.ticker.FuncFormatter(lambda x, pos: human_size(x * 1000, False)))
  return f


def do_accuracy_plot(frames, xid):
  f, axes = plt.subplots(len(frames) // 2, len(frames) // 2)
  axes = axes.flatten()
  for i, frame_id in enumerate(frames):
    sns.distplot(frames[frame_id][xid], ax=axes[i], label=frame_id)
    axes[i].annotate(frame_id, xy=(0.1, 0.85), size=10, textcoords='axes fraction')

  set_up_figure(f, names=(cols[xid], None), legend=False)
  f.suptitle('F-Score distribution')
  return f


def do_comparison_plot(frames, aid, bid):
  f, axes = plt.subplots(len(frames) // 2, len(frames) // 2)
  axes = axes.flatten()
  for i, frame_id in enumerate(frames):
    frame = frames[frame_id]
    sns.kdeplot(frame[aid], frame[bid], shade=True, ax=axes[i], label=frame_id)
    axes[i].annotate(frame_id, xy=(0.1, 0.85), size=10, textcoords='axes fraction')
  set_up_figure(f, names=(cols[aid], cols[bid]))
  return f


def latex_table(rows):
  result = ''
  for row in rows:
    row = np.array(row).tolist()
    result += ' & '.join([str(i) for i in row]) + ' \\\\\n'
  return result


def do_mfe_table(frames, xid):
  data = [i[xid] for i in frames.values()]
  data = pd.concat(data, axis=1)
  num_unique = data.apply(pd.Series.nunique, axis=1)
  num_unique = num_unique.value_counts().sort_index()
  num_unique /= num_unique.sum()
  table = [
    num_unique.index.tolist(),
    ['%.2f\\%%' % (i * 100) for i in num_unique.tolist()]
  ]
  print(latex_table(table))


def do_accuracy_table(frames, ids):
  table = [['Package'] + [
    '%s %s' % (cols[i], t) for i in ids for t in ['mean', 'SD']]]
  for frame_id, frame in frames.items():
    table.append([frame_id] + [
      '%.3f' % val for i in ids for val in (frame[i].mean(), frame[i].std())])
  print(latex_table(table))


def do_ttest(frames, xid):
  frame_keys = sorted(frames.keys())
  table = [['Package'] + frame_keys]
  for i, a_fid in enumerate(frame_keys):
    table.append([a_fid])
    for j, b_fid in enumerate(frame_keys):
      if i == j:
        table[-1].append('---')
        continue
      frameA, frameB = frames[a_fid][xid], frames[b_fid][xid]
      if i < j:
        frameA, frameB = frameB, frameA
      #t, p = scipy.stats.wilcoxon(frameA, frameB)
      t, p = scipy.stats.ttest_rel(frameA, frameB)
      entry = '%.3f' % t
      if i < j:
        entry = '%.3f' % p
      table[-1].append(entry)
  print(latex_table(table))


def read_frames(data):
  names = ['length', 'real', 'usersys', 'maxrss', 'fscore', 'ppv', 'sensitivity', 'mfe']
  return {
    name: pd.read_csv(
      filename, delimiter=' ', header=None,
      names=names)
    for name, filename in data.items()
    }


def savefig_local(dataset, name, f):
  save_figure(f, './benchmark_figures/%s_%s.png' % (dataset, name))


def do_perf_results(dataset_name, data):
  frames = read_frames(data)
  frames_by_length = {name: frame.groupby('length') for name, frame in frames.items()}

  # savefig_local(
  #   dataset_name, 'real',
  #   do_quantity_plot(frames_by_length, 'length', 'real'))
  savefig_local(
    dataset_name, 'maxrss',
    do_quantity_plot(frames_by_length, 'length', 'maxrss'))
  # f1, f2 = do_quantity_log_plot(frames_by_length, 'length', 'real')
  # savefig_local(dataset_name, 'real_loglog_scatter', f1)
  # savefig_local(dataset_name, 'real_loglog_bestfit', f2)

  f1, f2 = do_quantity_log_plot(frames_by_length, 'length', 'maxrss', logx=False)
  savefig_local(dataset_name, 'maxrss_loglog_scatter', f1)
  savefig_local(dataset_name, 'maxrss_loglog_bestfit', f2)


def do_accuracy_results(dataset_name, data):
  frames = read_frames(data)

  # savefig_local(dataset_name, 'fscore', do_accuracy_plot(frames, 'fscore'))
  # savefig_local(dataset_name, 'comparison', do_comparison_plot(frames, 'ppv', 'sensitivity'))
  # do_mfe_table(frames, 'mfe')
  # do_accuracy_table(frames, ['fscore', 'ppv', 'sensitivity'])
  # do_ttest(frames, 'fscore')


# ArchiveII
# do_accuracy_results('archiveii', {
#   'UNAFold': './benchmark_results/UNAFold_archiveii.results',
#   'RNAstructure': './benchmark_results/RNAstructure_archiveii.results',
#   'ViennaRNA-d3': './benchmark_results/ViennaRNAd3_archiveii.results',
#   'ViennaRNA-d2': './benchmark_results/ViennaRNAd2_archiveii.results'
# })

# Random
do_perf_results('random', {
  'UNAFold': './benchmark_results/UNAFold_random.results',
  'RNAstructure': './benchmark_results/RNAstructure_random.results',
  'ViennaRNA-d3': './benchmark_results/ViennaRNAd3_random.results',
  'ViennaRNA-d2': './benchmark_results/ViennaRNAd2_random.results'
})
