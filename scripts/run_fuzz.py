#!/usr/bin/env python3
import random

import math
import rna
import run
from common import read_file
from scrape import MemeVault

def gen_seqs(N):
  seqs = []
  for i in range(1 << (2 * N)):
    seq = ''
    while len(seq) != N:
      seq += 'GUAC'[i % 4]
      i //= 4
    seqs.append(seq)
  return seqs

def gen_rand_rna(N):
  return ''.join(['GUAC'[random.randint(0, 3)] for i in range(N)])

def gen_exhaustive():
  seqs = []
  for i in range(1, 10):
    seqs += gen_seqs(i)
  rnas = []
  for seq in seqs:
    rnas.extend(rna.generate_all_foldings(rna.RNA('meme', seq, [-1] * len(seq))))
  return rnas

def gen_random(N):
  seqs = []
  N2 = int(math.sqrt(N))
  for i in range(N2):
    seqs.append(gen_rand_rna(30))
  rnas = []
  for seq in seqs:
    rnas.extend(rna.generate_random_foldings(rna.RNA('meme', seq, [-1] * len(seq)), N2))
  return rnas

def read_fuzz_rnas(filename):
  return rna.rnas_from_db_list(read_file(filename))

def main():
  memevault = MemeVault('random')
  rnas = [r for r in memevault]
  print('Fuzzing with %d RNAs' % len(rnas))

  memerna = run.MemeRNA()
  rnastructure = run.RNAstructure()
  memerna_rnas = memerna.batch_efn(rnas)
  rnastructure_rnas = rnastructure.batch_efn(rnas)
  for i in range(len(rnas)):
    if memerna_rnas[0][i] != rnastructure_rnas[0][i]:
      print('Error on rna:\n%s\n %f - %f' % (rnas[i], memerna_rnas[0][i], rnastructure_rnas[0][i]))


if __name__ == "__main__":
  main()
