#!/usr/bin/env python3
import argparse
import os
import re

import sys

parser = argparse.ArgumentParser()
parser.add_argument(
  '-t', '--type', choices=['debug', 'asan', 'msan', 'ubsan', 'release', 'relwithdebinfo'],
  default='debug', required=False)
parser.add_argument('-c', '--use_clang', action='store_true', default=False, required=False)
parser.add_argument('-r', '--regenerate', action='store_true', default=False, required=False)
parser.add_argument('targets', nargs='*', type=str)
args = parser.parse_args()

defs = {
  'CMAKE_C_COMPILER': 'cc',
  'CMAKE_CXX_COMPILER': 'c++',
  'CMAKE_BUILD_TYPE': args.type
}
if args.use_clang:
  defs['CMAKE_C_COMPILER'] = 'clang'
  defs['CMAKE_CXX_COMPILER'] = 'clang++'

def run_command(cmd):
  res = os.system(cmd)
  if res != 0:
    sys.exit(1)


def should_regenerate_cmake():
  if not os.path.exists('build/CMakeCache.txt'):
    return True
  with open('build/CMakeCache.txt', 'r') as f:
    data = f.read()
    match = re.search('^CMAKE_BUILD_TYPE:STRING=(.*)$', data, re.M)
    if match is None or match.group(1) != defs['CMAKE_BUILD_TYPE']:
      return True
    match = re.search('^CMAKE_CXX_COMPILER:UNINITIALIZED=(.*)$', data, re.M)
    if match is None or match.group(1) != defs['CMAKE_CXX_COMPILER']:
      return True
  return False


regenerate = False
if not os.path.exists('build'):
  os.mkdir('build')
  regenerate = True
regenerate = regenerate or should_regenerate_cmake() or args.regenerate

os.chdir('build')
if regenerate:
  print('Regenerating cmake files.')
  def_str = ' '.join('-D %s=\'%s\'' % (i, k) for i, k in defs.items())
  run_command('cmake %s ../' % def_str)
run_command('make -j16 %s' % ' '.join(args.targets))
os.chdir('../')
