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
  f.suptitle('%s vs %s' % tuple(reversed(names)))
  for ax in f.get_axes():
    set_up_axis(ax, names, legend)


def save_figure(f, name):
  f.tight_layout()
  f.set_size_inches(12, 8)
  f.savefig(name, dpi=100)


def latex_table(rows):
  result = ''
  for row in rows:
    row = np.array(row).tolist()
    result += ' & '.join([str(i) for i in row]) + ' \\\\\n'
  return result


def savefig_local(dataset, name, f):
  save_figure(f, './build/benchmark_figures/%s_%s.png' % (dataset, name))


def get_best_factors(n):
  best = 1
  for i in range(1, int(np.floor(np.sqrt(n))) + 1):
    if n % i == 0:
      best = i
  return best, n // best


def get_best_splitting(n):
  a = get_best_factors(n)
  if a[0] == 1 and n != 1:
    return get_best_factors(n + 1)
  return a


def get_subplot_grid(n, sharex=False, sharey=False):
  factors = get_best_splitting(n)
  if factors == (1, 1):
    f, axes = plt.subplots(n, sharey=sharey, sharex=sharex)
  else:
    f, axes = plt.subplots(factors[0], factors[1], sharey=sharey, sharex=sharex)
    axes = axes.flatten()
  if n == 1:
    axes = [axes]
  return f, axes
