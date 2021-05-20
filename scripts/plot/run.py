#!/usr/bin/env python3
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
import os
import sys

sys.path.append('scripts')

import pandas as pd
from scripts.plot.plot_common import savefig_local
from scripts.plot.subopt_performance import subopt_distribution
from scripts.plot.fold_performance import fold_perf_results
from scripts.plot.fold_accuracy import fold_accuracy_results
from scripts.plot.load_data import read_fold_dataset, read_subopt_dataset, load_subset_file

import seaborn as sns

PREFIX = '../results'


def generate_filename_map(dataset_name, enable):
  d = {
    'ViennaRNA-d2': os.path.join(PREFIX, 'ViennaRNA-d2_%s.results' % dataset_name),
    'ViennaRNA-d2-sorted': os.path.join(PREFIX, 'ViennaRNA-d2-sorted_%s.results' % dataset_name),
    'SparseMFEFold': os.path.join(PREFIX, 'SparseMFEFold_%s.results' % dataset_name),
    'memerna': os.path.join(PREFIX, 'memerna_%s.results' % dataset_name),
    'memerna-sorted': os.path.join(PREFIX, 'memerna-sorted_%s.results' % dataset_name),
    'RNAstructure': os.path.join(PREFIX, 'RNAstructureDistribution_%s.results' % dataset_name),
    'RNAstructure-mod': os.path.join(PREFIX, 'RNAstructureHarness_%s.results' % dataset_name),
    'ViennaRNA-d3': os.path.join(PREFIX, 'ViennaRNA-d3_%s.results' % dataset_name),
    'ViennaRNA-d3-sorted': os.path.join(PREFIX, 'ViennaRNA-d3-sorted_%s.results' % dataset_name),
    'SJSViennaMPI': os.path.join(PREFIX, 'SJSViennaMPI_%s.results' % dataset_name),
    'SJSViennaMPI-sorted': os.path.join(PREFIX, 'SJSViennaMPI-sorted_%s.results' % dataset_name),
    'UNAFold': os.path.join(PREFIX, 'UNAFold_%s.results' % dataset_name)
  }
  return {k: d[k] for k in enable}


FAST_FOLDERS = ['ViennaRNA-d2', 'SparseMFEFold', 'memerna', 'ViennaRNA-d3']
PERF_FOLDERS = FAST_FOLDERS + ['RNAstructure', 'UNAFold']
ACCURACY_FOLDERS = PERF_FOLDERS + ['RNAstructure-mod']


def fold_accuracy_graphs():
  fold_archiveii_domains_ds = read_fold_dataset(
    'archiveii_domains',
    generate_filename_map('archiveii', ACCURACY_FOLDERS),
    os.path.join(PREFIX, 'archiveii_domains_dedup.subset')
  )
  fold_accuracy_results(fold_archiveii_domains_ds)


def fold_graphs():
  fold_random_all_ds = read_fold_dataset(
    'random',
    generate_filename_map('random', PERF_FOLDERS),
    os.path.join(PREFIX, 'random_all.subset')
  )
  fold_random_all_fast_ds = read_fold_dataset(
    'random_fast',
    generate_filename_map('random', FAST_FOLDERS),
    os.path.join(PREFIX, 'random_all.subset')
  )
  fold_random_large_all_ds = read_fold_dataset(
    'random_large',
    generate_filename_map('random_large', FAST_FOLDERS + ['UNAFold']),
    os.path.join(PREFIX, 'random_large_all.subset')
  )
  print('Loaded data')
  fold_perf_results(fold_random_all_ds, True, True)
  fold_perf_results(fold_random_all_fast_ds, True, False)
  fold_perf_results(fold_random_large_all_ds, True, True)


DELTAS = [1, 2, 3, 4, 5, 6, 10, 11, 12, 13]


def subopt_graphs():
  ALL_SUBOPTS = ['ViennaRNA-d2', 'RNAstructure', 'ViennaRNA-d3', 'ViennaRNA-d2-sorted',
                 'ViennaRNA-d3-sorted', 'SJSViennaMPI', 'SJSViennaMPI-sorted', 'memerna',
                 'memerna-sorted']
  # ALL_SUBOPTS = ['RNAstructure',  'memerna', 'SJSViennaMPI', 'SJSViennaMPI-sorted']
  # ALL_SUBOPTS = ['ViennaRNA-d2', 'ViennaRNA-d2-sorted', 'ViennaRNA-d3', 'ViennaRNA-d3-sorted',]
  # ALL_SUBOPTS = ['ViennaRNA-d2', 'RNAstructure', 'ViennaRNA-d3', 'SJSViennaMPI', 'memerna']
  # ALL_SUBOPTS = ['ViennaRNA-d2', 'ViennaRNA-d2-sorted', 'memerna-sorted', 'memerna']
  # ALL_SUBOPTS = ['memerna', 'ViennaRNA-d2']
  deltas = [1, 2, 3, 4, 5, 6, 10, 11, 12, 13]

  all_ds = []
  for delta in deltas:
    subopt_random_all_ds = read_subopt_dataset(
      'random_subopt_%d' % delta,
      generate_filename_map('random_subopt_%d' % delta, ALL_SUBOPTS),
      os.path.join(PREFIX, 'random_all.subset')
    )
    all_ds.append(subopt_random_all_ds.fmap)
    # subopt_perf_results(subopt_random_all_ds)
  savefig_local(
    'subopt_distribution', 'numstruc',
    subopt_distribution(all_ds, ALL_SUBOPTS))
  # subopt_max_ds = read_subopt_dataset(
  #   'random_subopt_max',
  #   generate_filename_map('random_subopt_max', ['memerna', 'ViennaRNA-d2']),
  #   os.path.join(PREFIX, 'random_all.subset')
  # )
  # subopt_perf_results(subopt_max_ds)


def extract_highest_structure_count_rows(subopt_name):
  all_ds = []
  cols = ['name', 'run', 'length', 'real', 'usersys', 'maxrss', 'numstruc']
  for delta in DELTAS:
    fmap = generate_filename_map('random_subopt_%d' % delta, [subopt_name])
    frame = pd.read_csv(fmap[subopt_name], delimiter=' ', header=None, names=cols)
    all_ds.append(frame.groupby(frame['name']))
  subset = load_subset_file(os.path.join(PREFIX, 'random_all.subset'))

  def lookup(ds):
    try:
      return ds[1].get_group(s)['numstruc'].iloc[0]
    except:
      return -1

  all_rows = []
  for s in subset:
    idx, val = max(enumerate(all_ds), key=lookup)
    all_rows += val.get_group(s).values.tolist()
  max_frame = pd.DataFrame(all_rows, columns=cols)
  max_frame.to_csv('test.results', sep=' ', header=False, index=False)


def main():
  sns.set(color_codes=True)
  # fold_accuracy_graphs()
  # fold_graphs()
  # subopt_graphs()
  # extract_highest_structure_count_rows('ViennaRNA-d2-sorted')


if __name__ == '__main__':
  main()
