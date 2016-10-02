from plot.plot_common import savefig_local, do_quantity_plot, do_quantity_log_plot


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
