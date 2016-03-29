#!/usr/bin/env python3
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--type', choices=['debug', 'asan', 'msan', 'ubsan'], default='debug', required=False)
parser.add_argument('-c', '--use_clang', action='store_true', default=False, required=False)
args = parser.parse_args()

cc = 'cc'
cxx = 'c++'
if args.use_clang:
  cc = 'clang'
  cxx = 'clang++'

def run_command(cmd):
  print(cmd)
  os.system(cmd)

if not os.path.exists('build'):
  os.mkdir('build')
os.chdir('build')
run_command('cmake -D CC=%s -D CXX=%s -D CMAKE_BUILD_TYPE=%s ../' % (cc, cxx, args.type))
run_command('make -j16')
os.chdir('../')
