#!/usr/bin/env python3
# Copyright 2016 E.
import argparse
import os
import sys
import shutil
from pathlib import Path

parser = argparse.ArgumentParser()
# Build environment options
parser.add_argument('-p', '--prefix', type=str, default=os.path.join(Path.home(), "bin"),
  required=False, help='Where to place build directory')
parser.add_argument(
  '-t', '--type', choices=['debug', 'asan', 'ubsan', 'release', 'relwithdebinfo'],
  default='debug', required=False)
parser.add_argument('-c', '--use-clang', action='store_true', default=False, required=False)
parser.add_argument('-a', '--use-afl', action='store_true', default=False, required=False)
parser.add_argument('-r', '--regenerate', action='store_true', default=False, required=False)
parser.add_argument('--compilers', type=str, nargs=2, required=False)
parser.add_argument('targets', nargs='*', type=str)

# Memerna configuration options:
parser.add_argument('-m', '--use-mpfr', action='store_true', default=False, required=False)

# Misc options:
parser.add_argument('-d', '--dry', action='store_true', default=False, required=False)

args = parser.parse_args()

def run_command(cmd):
  if args.dry:
    print(cmd)
  else:
    res = os.system(cmd)
    if res != 0:
      sys.exit(1)

if bool(args.use_clang) + bool(args.use_afl) + bool(args.compilers) > 1:
  parser.error('At most one compiler related flag allowed simultaneously')

compilers = ('cc', 'c++')
env = [
  # Add stack protector etc to catch non-crashing memory bugs.
  'AFL_HARDEN=1',
  # Find more paths
  'AFL_LLVM_LAF_ALL=1'
]

if args.use_clang:
  compilers = ('clang', 'clang++')
elif args.use_afl:
  compilers = ('afl-clang-lto', 'afl-clang-lto++')
elif args.compilers:
  compilers = args.compilers

defs = {
  'CMAKE_C_COMPILER': compilers[0],
  'CMAKE_CXX_COMPILER': compilers[1],
  'CMAKE_BUILD_TYPE': args.type,
  'PARTITION_MPFR': 'OFF'
}

if args.use_mpfr:
  defs['PARTITION_MPFR'] = 'ON'

build_dir = os.path.join(args.prefix, 'memerna', defs['CMAKE_CXX_COMPILER'] + '-' + defs['CMAKE_BUILD_TYPE'])
if defs['PARTITION_MPFR'] == 'ON':
  build_dir += '-mpfr'
regenerate = args.regenerate

if regenerate and os.path.exists(build_dir):
  shutil.rmtree(build_dir)

print(build_dir)
if not os.path.exists(build_dir):
  os.makedirs(build_dir)
  regenerate = True

proj_dir = os.getcwd()

os.chdir(build_dir)
if regenerate:
  print('Regenerating cmake files.')
  def_str = ' '.join('-D %s=\'%s\'' % (i, k) for i, k in defs.items())
  run_command('cmake %s %s' % (def_str, proj_dir))
run_command('%s make -j$(nproc) %s' % (' '.join(env), ' '.join(args.targets)))
os.chdir('../../')
