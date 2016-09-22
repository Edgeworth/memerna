import os
import resource
import subprocess
import sys


def float_fmt(f):
  return ('%.2f' % f).rstrip('0').rstrip('.')


def human_size(b, binary=True):
  units = ['B', 'KiB', 'MiB', 'GiB']
  base = 1024
  if not binary:
    units = ['B', 'KB', 'MB', 'GB']
    base = 1000
  for unit in units[:-1]:
    if abs(b) < base:
      return '%s %s' % (float_fmt(b), unit)
    b /= base
  return '%s %s' % (float_fmt(b), units[-1])


class ProcessResults:
  def __init__(self, stdout, stderr, ret, real, usersys, maxrss):
    self.maxrss = maxrss
    self.usersys = usersys
    self.real = real
    self.ret = ret
    self.stdout = stdout
    self.stderr = stderr

  def __str__(self):
    return '%.2fs, %s ' % (self.real, human_size(self.maxrss))


def try_command(*cmd, record_stdout=False, input=None, limits=None):
  if isinstance(record_stdout, str):
    stdout = open(record_stdout, 'w')
  else:
    stdout = subprocess.PIPE if record_stdout else subprocess.DEVNULL
  stdin = subprocess.PIPE if input else None
  if input:
    input = input.encode('utf-8')
  cmd = ['/usr/bin/time', '-f', '%e %U %S %M'] + list(cmd)

  def pre_exec():
    if limits:
      if limits[0]:
        resource.setrlimit(resource.RLIMIT_AS, (limits[0] * 1024, limits[0] * 1024))
      if limits[1]:
        resource.setrlimit(resource.RLIMIT_CPU, (limits[1], limits[1]))

  with subprocess.Popen(
      cmd, shell=False, stdin=stdin, stdout=stdout,
      stderr=subprocess.PIPE, preexec_fn=pre_exec) as proc:
    stdout_data, stderr_data = proc.communicate(input=input)
    ret = proc.wait()
    last_line = stderr_data.decode('UTF-8').strip().split('\n')[-1].split(' ')
    real, user, sys, maxrss = [float(i) for i in last_line]
    if stdout_data:
      stdout_data = stdout_data.decode('UTF-8')
    if isinstance(record_stdout, str):
      stdout.flush()
      stdout.close()
    return ProcessResults(stdout_data, stderr_data, ret, real, user + sys, maxrss * 1024)


def run_command(*cmd, record_stdout=False, input=None, limits=None):
  res = try_command(*cmd, record_stdout=record_stdout, input=input, limits=limits)
  if res.ret:
    print('Running `%s\' failed with ret code %d.\nStderr:\n%s\n' % (
      cmd, res.ret, res.stderr.decode('utf-8')))
    sys.exit(1)
  return res


def fix_path(path):
  return os.path.abspath(os.path.expanduser(path))


def read_file(name):
  return open(name, encoding='utf-8').read()


def write_file(name, data):
  with open(name, 'w') as f:
    f.write(data)
