#!/usr/bin/env python3
import argparse
import os
import tempfile

from common import run_command
from convert import vienna_to_ct

def run_rnastructure(rnastrucuture_dir, seq, pairs):
  print("Running RNAstructure")
  with tempfile.NamedTemporaryFile() as f, tempfile.NamedTemporaryFile() as out:
    f.write(vienna_to_ct(seq, pairs).encode('UTF-8'))
    f.flush()
    os.putenv('DATAPATH', '%s/data_tables' % rnastrucuture_dir)
    run_command('%s/exe/efn2' % rnastrucuture_dir, '-s', f.name, out.name)
    print('RNAstructure: %s' % out.read().strip().decode('UTF-8'))

def run_memerna(memerna_loc, seq, pairs):
  print("Running memerna")
  run_command(memerna_loc, seq, pairs)

def run_mrna(mrna_loc, seq, pairs):
  print("Running mrna")


def run_vienna_rna(vienna_rna, seq, pairs):
  print("Running ViennaRNA")


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('filename')
  parser.add_argument('--rnastructure-dir')
  parser.add_argument('--memerna-loc')
  parser.add_argument('--mrna-loc')
  parser.add_argument('--vienna-rna-loc')
  args = parser.parse_args()

  seq, pairs = open(args.filename).readlines()
  seq = seq.strip()
  pairs = pairs.strip()

  if args.rnastructure_dir:
    run_rnastructure(args.rnastructure_dir, seq, pairs)
  if args.memerna_loc:
    run_memerna(args.memerna_loc, seq, pairs)
  # run_mrna(args.mrna_loc, seq, pairs)
  # run_vienna_rna(args.vienna_rna, seq, pairs)

if __name__ == "__main__":
  main()
