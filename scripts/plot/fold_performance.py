# Copyright 2016 Eliot Courtney.
from scripts.plot.plot_common import savefig_local, do_quantity_plot, do_quantity_log_plot


def fold_perf_results(ds, do_normal_plots=True, do_log_plots=True):
  fmap_by_len = {name: frame.groupby('length') for name, frame in ds.fmap.items()}

  if do_normal_plots:
    savefig_local(
      ds.name, 'real',
      do_quantity_plot(fmap_by_len, 'length', 'real'))
    savefig_local(
      ds.name, 'maxrss',
      do_quantity_plot(fmap_by_len, 'length', 'maxrss'))
  if do_log_plots:
    f1 = do_quantity_log_plot(fmap_by_len, 'length', 'real')
    savefig_local(ds.name, 'real_loglog_scatter', f1)

    f1 = do_quantity_log_plot(fmap_by_len, 'length', 'maxrss')
    savefig_local(ds.name, 'maxrss_loglog_scatter', f1)
