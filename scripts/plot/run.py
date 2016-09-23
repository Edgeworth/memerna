#!/usr/bin/env python3
import os
import sys

sys.path.append('scripts')

from plot.fold_performance import fold_perf_results
from plot.fold_accuracy import fold_accuracy_results
from plot.load_data import read_fold_dataset

import seaborn as sns

PREFIX = '../results'


def generate_filename_map(dataset_name, rnastructure=True, vd3=True, unafold=True):
  d = {
    'ViennaRNA-d2': os.path.join(PREFIX, 'ViennaRNA-d2_%s.results' % dataset_name),
    'SparseMFEFold': os.path.join(PREFIX, 'SparseMFEFold_%s.results' % dataset_name),
    'memerna': os.path.join(PREFIX, 'MemeRNA_%s.results' % dataset_name)
  }
  if rnastructure:
    d['RNAstructure'] = os.path.join(PREFIX, 'RNAstructureDistribution_%s.results' % dataset_name)
    d['RNAstructure-mod'] = os.path.join(PREFIX, 'RNAstructureHarness_%s.results' % dataset_name)
  if vd3:
    d['ViennaRNA-d3'] = os.path.join(PREFIX, 'ViennaRNA-d3_%s.results' % dataset_name)
  if unafold:
    d['UNAFold'] = os.path.join(PREFIX, 'UNAFold_%s.results' % dataset_name)
  return d


def main():
  sns.set(color_codes=True)
  # fold_test_ds = read_fold_dataset('test', {'test_prog': 'test.results'}, 'test.subset')
  # fold_archiveii_all_ds = read_fold_dataset(
  #   'archiveii',
  #   generate_filename_map('archiveii'),
  #   os.path.join(PREFIX, 'archiveii_all_dedup.subset')
  # )
  fold_archiveii_domains_ds = read_fold_dataset(
    'archiveii_domains',
    generate_filename_map('archiveii'),
    os.path.join(PREFIX, 'archiveii_domains_dedup.subset')
  )
  fold_random_all_ds = read_fold_dataset(
    'random',
    generate_filename_map('random'),
    os.path.join(PREFIX, 'random_all.subset')
  )
  fold_random_all_fast_ds = read_fold_dataset(
    'random_fast',
    generate_filename_map('random', False, False, False),
    os.path.join(PREFIX, 'random_all.subset')
  )
  fold_random_large_all_ds = read_fold_dataset(
    'random_large',
    generate_filename_map('random_large', False),
    os.path.join(PREFIX, 'random_large_all.subset')
  )
  print('Loaded data')
  # fold_perf_results(fold_test_ds)
  fold_perf_results(fold_random_all_ds)
  fold_perf_results(fold_random_all_fast_ds)
  fold_perf_results(fold_random_large_all_ds)
  fold_accuracy_results(fold_archiveii_domains_ds)


if __name__ == '__main__':
  main()
