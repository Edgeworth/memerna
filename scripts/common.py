import os
import subprocess

import numpy as np


class BenchmarkResults:
  def __init__(self, real, user, sys, maxrss):
    self.real = real
    self.user = user
    self.sys = sys
    self.maxrss = maxrss

def run_command(*args, input=None):
  print("Running `%s'" % ' '.join(args))
  return subprocess.run(
    args, input=input, stdout=subprocess.PIPE,
    stderr=subprocess.PIPE, check=True)


def benchmark_command(*args, num_runs=10):
  cum_results = np.zeros(4)
  stdouts = []
  for i in range(num_runs):
    res = run_command('/usr/bin/time', '-f', '%e %U %S %M', *args)
    stdouts.append(res.stdout.decode('UTF-8'))
    t = res.stderr.decode('UTF-8')
    cum_results += np.array([float(i) for i in t.split(' ')])
  cum_results /= num_runs
  return (
    BenchmarkResults(
      real=cum_results[0], user=cum_results[1],
      sys=cum_results[2], maxrss=cum_results[3]),
    stdouts
  )


def fix_path(path):
  return os.path.abspath(os.path.expanduser(path))


def read_file(name):
  return open(name, encoding='utf-8').read()


def write_file(name, data):
  with open(name, 'w') as f:
    f.write(data)
