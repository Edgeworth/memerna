import os
import shlex
import sys


def run_command(cmd, *args):
  cmd = '%s %s' % (shlex.quote(cmd), ' '.join(shlex.quote(i) for i in args))
  print("Running `%s'" % cmd)
  if os.system(cmd) != 0:
    sys.exit(1)


def read_file(name):
  return open('orig_data/%s' % name, encoding='utf-8').read()


def write_file(name, data):
  with open('data/%s' % name, 'w') as f:
    f.write(data)
