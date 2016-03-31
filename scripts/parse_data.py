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
      i += 5
      idx += 1
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
