import numpy as np

from matplotlib import pyplot as plt

EP = 1e-2


def set_up_axis(ax, names, legend):
  for i, axis in [(0, ax.xaxis), (1, ax.yaxis)]:
    if names[i] is not None:
      axis.set_label_text('%s' % names[i])
  if legend:
    ax.legend(loc=0, fontsize='medium')


def set_up_figure(f, names=None, legend=True):
  f.suptitle('%s vs %s' % tuple(reversed(names)), y=1.00)
  for ax in f.get_axes():
    set_up_axis(ax, names, legend)


def save_figure(f, name):
  f.tight_layout()
  f.savefig(name, dpi=150)


def savefig_local(dataset, name, f):
  save_figure(f, './build/benchmark_figures/%s_%s.png' % (dataset, name))
  plt.close(f)


def latex_table(rows):
  result = ''
  for row in rows:
    row = np.array(row).tolist()
    result += ' & '.join([str(i) for i in row]) + ' \\\\\n'
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
