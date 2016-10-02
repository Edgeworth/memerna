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

from plot.subopt_performance import subopt_perf_results
from plot.fold_performance import fold_perf_results
from plot.fold_accuracy import fold_accuracy_results
from plot.load_data import read_fold_dataset, read_subopt_dataset

import seaborn as sns

PREFIX = '../results'


def generate_filename_map(dataset_name, enable):
  # TODO SJSVienna, SJSViennaMPI
  d = {
    'ViennaRNA-d2': os.path.join(PREFIX, 'ViennaRNA-d2_%s.results' % dataset_name),
    'ViennaRNA-d2-sorted': os.path.join(PREFIX, 'ViennaRNA-d2-sorted_%s.results' % dataset_name),
    'SparseMFEFold': os.path.join(PREFIX, 'SparseMFEFold_%s.results' % dataset_name),
    'memerna': os.path.join(PREFIX, 'MemeRNA_%s.results' % dataset_name),
    'RNAstructure': os.path.join(PREFIX, 'RNAstructureDistribution_%s.results' % dataset_name),
    'RNAstructure-mod': os.path.join(PREFIX, 'RNAstructureHarness_%s.results' % dataset_name),
    'ViennaRNA-d3': os.path.join(PREFIX, 'ViennaRNA-d3_%s.results' % dataset_name),
    'ViennaRNA-d3-sorted': os.path.join(PREFIX, 'ViennaRNA-d3-sorted_%s.results' % dataset_name),
    'SJSViennaMPI': os.path.join(PREFIX, 'SJSViennaMPI_%s.results' % dataset_name),
    'SJSViennaMPI-sorted': os.path.join(PREFIX, 'SJSViennaMPI-sorted_%s.results' % dataset_name),
    'UNAFold': os.path.join(PREFIX, 'UNAFold_%s.results' % dataset_name)
  }
  return {k:d[k] for k in enable}


def fold_graphs():
  ALL_FOLDERS = ['ViennaRNA-d2', 'SparseMFEFold', 'memerna',
                 'RNAstructure', 'RNAstructure-mod', 'ViennaRNA-d3', 'UNAFold']
  FAST_FOLDERS = ['ViennaRNA-d2', 'SparseMFEFold', 'memerna', 'ViennaRNA-d3']
  fold_archiveii_domains_ds = read_fold_dataset(
    'archiveii_domains',
    generate_filename_map('archiveii', ALL_FOLDERS),
    os.path.join(PREFIX, 'archiveii_domains_dedup.subset')
  )
  fold_random_all_ds = read_fold_dataset(
    'random',
    generate_filename_map('random', ALL_FOLDERS),
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
  fold_perf_results(fold_random_all_ds)
  fold_perf_results(fold_random_all_fast_ds)
  fold_perf_results(fold_random_large_all_ds)
  fold_accuracy_results(fold_archiveii_domains_ds)


def main():
  sns.set(color_codes=True)
  # TODO SJS Vienna
  ALL_SUBOPTS = ['ViennaRNA-d2', 'RNAstructure', 'RNAstructure-mod', 'ViennaRNA-d3',
                 'ViennaRNA-d2-sorted', 'ViennaRNA-d3-sorted']
  # ALL_SUBOPTS = ['memerna']
  deltas = [1, 2, 3, 4, 5, 6, 10, 11, 12, 13]
  for delta in deltas:
    subopt_random_all_ds = read_subopt_dataset(
      'random_subopt_%d' % delta,
      generate_filename_map('random_subopt_%d' % delta, ALL_SUBOPTS),
      os.path.join(PREFIX, 'random_all.subset')
    )
    subopt_perf_results(subopt_random_all_ds)


if __name__ == '__main__':
  main()
