#!/usr/bin/env python3
import argparse
import os

import sys

parser = argparse.ArgumentParser()
parser.add_argument(
  '-t', '--type', choices=['debug', 'asan', 'msan', 'ubsan', 'release', 'relwithdebinfo'],
  default='debug', required=False)
parser.add_argument('-c', '--use_clang', action='store_true', default=False, required=False)
parser.add_argument('targets', nargs='*', type=str)
args = parser.parse_args()

cc = 'cc'
cxx = 'c++'
defs = {'CMAKE_C_COMPILER': 'cc', 'CMAKE_CXX_COMPILER': 'c++'}
if args.use_clang:
  defs['CMAKE_C_COMPILER'] = 'clang'
  defs['CMAKE_CXX_COMPILER'] = 'clang++'
  defs['CMAKE_CXX_FLAGS_MSAN'] = '-fsanitize-memory-track-origins  -fsanitize-memory-use-after-dtor'
  defs['CMAKE_CXX_FLAGS_UBSAN'] = '-fsanitize=unsigned-integer-overflow'

def run_command(cmd):
  res = os.system(cmd)
  if res != 0:
    sys.exit(1)


if not os.path.exists('build'):
  os.mkdir('build')
os.chdir('build')
def_str = ' '.join('-D %s=\'%s\'' % (i, k) for i, k in defs.items())
run_command('cmake %s -D CMAKE_BUILD_TYPE=%s ../' % (def_str, args.type))
run_command('make -j16 %s' % ' '.join(args.targets))
os.chdir('../')
