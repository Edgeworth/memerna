# Copyright 2016, E.
#
# This file is part of memerna.
#
# memerna is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with memerna.
# If not, see <http://www.gnu.org/licenses/>.
from plot.plot_common import savefig_local, do_quantity_plot, do_quantity_log_plot


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
