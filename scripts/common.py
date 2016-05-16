import logging
import os
import subprocess

import numpy as np

log = logging.getLogger('')


def human_size(b):
  for unit in ['B', 'KiB', 'MiB']:
    if abs(b) < 1024:
      return '%.2f %s' % (b, unit)
    b /= 1024
  return "%.2f %s" % (b, 'GiB')


class BenchmarkResults:
  def __init__(self, real, usersys, maxrss):
    self.real = real
    self.usersys = usersys
    self.maxrss = maxrss

  def __str__(self):
    return '%.2f real, %.2f usersys %s max rss' % (
      self.real, self.usersys, human_size(self.maxrss))

  def to_np_array(self):
    return np.array([self.real, self.usersys, self.maxrss])

  @staticmethod
  def from_np_array(array):
    return BenchmarkResults(*array.tolist())

  @staticmethod
  def combine_benchmarks(benchmarks):
    num_include = 0
    total = np.zeros(3)
    for i in benchmarks:
      if i.real < 0 or i.usersys < 0:
        continue
      total += i.to_np_array()
      num_include += 1
    if num_include < len(benchmarks):
      log.debug('Removed %d results with negative time' % (len(benchmarks) - num_include))
    return BenchmarkResults.from_np_array(total / num_include)


def run_command(*args, input=None):
  log.debug("Running `%s'" % ' '.join(args))
  return subprocess.run(
    args, input=input, stdout=subprocess.PIPE,
    stderr=subprocess.PIPE, check=True)


def benchmark_command(*args, input=None):
  stdouts = []
  res = run_command('/usr/bin/time', '-f', '%e %U %S %M', *args, input=input)
  real, user, sys, maxrss = [float(i) for i in res.stderr.decode('UTF-8').split(' ')]
  return (
    BenchmarkResults(real=real, usersys=user + sys, maxrss=maxrss),
    res.stdout.decode('UTF-8')
  )


def fix_path(path):
  return os.path.abspath(os.path.expanduser(path))


def read_file(name):
  return open(name, encoding='utf-8').read()


def write_file(name, data):
  with open(name, 'w') as f:
    f.write(data)
