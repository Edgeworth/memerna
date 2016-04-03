#!/usr/bin/env python3
import re

MAX = 1000000000
ORDER = 'ACGU'


def parse_number(val):
  if val == '.':
    return MAX
  else:
    assert(float(val) * 10 == int(val.replace('.', '')))
    return int(val.replace('.', ''))


def parse_2x2_file(data):
  lines = [i.strip() for i in re.sub(r' +', ' ', data).split('\n')]
  output = ''
  idx = 0
  for i in range(len(lines)):
    if re.match(r'(\s*3\' <-- 5\'\s*){4}', lines[i]):
      matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 1:i + 5]]
      for m in range(4):
        for r in range(4):
          for c in range(4):
            val = parse_number(matrix_lines[r][m * 4 + c])
            output += '%s%s%s%s %d\n' % (ORDER[idx], ORDER[r], ORDER[c], ORDER[m], val)
      idx += 1
  return output

def parse_1x1_internal_loop(data):
  lines = [i.strip() for i in re.sub(r' +', ' ', data).split('\n')]
  output = ''
  for i in range(len(lines)):
    if re.match(r'(\s*5\' --> 3\'\s*){6}', lines[i]):
      t3prime = re.findall(r'[GUAC]', lines[i+2])
      t5prime = re.findall(r'[GUAC]', lines[i+3])
      assert(len(t3prime) == 12)
      assert(len(t5prime) == 12)
      matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 6:i + 10]]
      for m in range(6):
        for r in range(4):
          for c in range(4):
            val = parse_number(matrix_lines[r][m * 4 + c])
            output += '%s%s%s%s%s%s %d\n' % (
                t3prime[2 * m], ORDER[r], t3prime[2 * m + 1], t5prime[2 * m + 1], ORDER[c], t5prime[2 * m], val)
  return output

def parse_1x2_internal_loop(data):
  lines = [i.strip() for i in re.sub(r' +', ' ', data).split('\n')]
  output = ''
  for i in range(len(lines)):
    if re.match(r'(\s*5\' --> 3\'\s*){6}', lines[i]):
      t3prime = re.findall(r'[GUAC]', lines[i+2])
      t5prime = re.findall(r'[GUAC]', lines[i+3])
      assert(len(t3prime) == 12)
      assert(len(t5prime) == 12)
      extra = re.findall(r'[GUAC]', lines[i+4])
      assert(len(extra) == 6)
      matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 6:i + 10]]
      for m in range(6):
        for r in range(4):
          for c in range(4):
            val = parse_number(matrix_lines[r][m * 4 + c])
            output += '%s%s%s%s%s%s%s %d\n' % (
              t3prime[2 * m], ORDER[r], t3prime[2 * m + 1], t5prime[2 * m + 1], extra[m], ORDER[c], t5prime[2 * m], val)
  return output

def parse_2x2_internal_loop(data):
  lines = [i.strip() for i in re.sub(r' +', ' ', data).split('\n')]
  output = ''
  # Skip first example.
  for i in range(15, len(lines)):
    if re.match(r'\s*5\' ------> 3\'\s*', lines[i]):
      t3prime = re.findall(r'[GUAC]', lines[i+1])
      t5prime = re.findall(r'[GUAC]', lines[i+2])
      assert(len(t3prime) == 2)
      assert(len(t5prime) == 2)
      matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 4:i + 20]]
      for x1 in range(4):
        for x2 in range(4):
          for y1 in range(4):
            for y2 in range(4):
              val = parse_number(matrix_lines[4 * x1 + x2][4 * y1 + y2])
              output += '%s%s%s%s%s%s%s%s %d\n' % (
                  t3prime[0], ORDER[x1], ORDER[y1], t3prime[1],
                  t5prime[1], ORDER[y2], ORDER[x2], t5prime[0], val)
  return output

def parse_map_file(data):
  m = re.findall(r'([GUAC]+)\s*(\S+)', data)
  return ''.join('%s %s\n' % (i[0], i[1].replace('.', '')) for i in m)

def parse_loop_file(data):
  m = re.findall(r'(\d+)\s+([0-9.\-+]+)\s+([0-9.\-+]+)\s+([0-9.\-+]+)', data)
  internal, bulge, hairpin = '', '', ''
  for i in m:
    internal += '%s %d\n' % (i[0], parse_number(i[1]))
    bulge += '%s %d\n' % (i[0], parse_number(i[2]))
    hairpin += '%s %d\n' % (i[0], parse_number(i[3]))
  return (internal, bulge, hairpin)

# Outputs AXYA number
def parse_stack_txt(data):
  return parse_2x2_file(data)


def parse_terminal_txt(data):
  return parse_2x2_file(data)


def read_file(name):
  return open('orig_data/%s' % name, encoding='utf-8').read()


def write_file(name, data):
  with open('data/%s' % name, 'w') as f:
    f.write(data)


write_file(
  'hairpin.data',
  parse_map_file(read_file('triloop.txt')) +
  parse_map_file(read_file('tloop.txt')) +
  parse_map_file(read_file('hexaloop.txt')))
write_file('stacking.data', parse_stack_txt(read_file('stack.txt')))
write_file('terminal.data', parse_terminal_txt(read_file('tstack.txt')))
internal, bulge, hairpin = parse_loop_file(read_file('loop.txt'))
write_file('internal_initiation.data', internal)
write_file('bulge_initiation.data', bulge)
write_file('hairpin_initiation.data', hairpin)
write_file('internal_1x1.data', parse_1x1_internal_loop(read_file('int11.txt')))
write_file('internal_1x2.data', parse_1x2_internal_loop(read_file('int21.txt')))
write_file('internal_2x2.data', parse_2x2_internal_loop(read_file('int22.txt')))
