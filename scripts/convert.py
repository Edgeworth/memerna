#!/usr/bin/env python3
import argparse
import re


def vienna_to_ct(seq, pairs):
  seq, pairs = seq.strip(), pairs.strip()
  ct = ['%d %s' % (len(seq), 'somerna')]

  m = {}
  s = []
  for i, v in enumerate(pairs):
    if v == '(':
      s.append(i)
    elif v == ')':
      p = s[-1]
      s.pop(-1)
      m[p] = i
      m[i] = p

  for i, v in enumerate(seq):
    ct.append('%d %s %d %d %d %d' % (i + 1, v, i, i + 2, m.get(i, -1) + 1, i + 1))

  return '\n'.join(ct)


def ct_to_vienna(ct):
  ct = re.sub(r' +', ' ', ct)

  ct = [i.split() for i in ct.split('\n')[1:] if i]
  print(ct)
  sequence = ''.join(i[1] for i in ct)
  pairs = ''
  for i, v in enumerate(ct):
    idx = int(v[4])
    if idx == 0:
      pairs += '.'
    elif idx < i:
      pairs += ')'
    else:
      pairs += '('
  return sequence, pairs


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('filename')
  args = parser.parse_args()

  seq, pairs = open(args.filename).readlines()
  print(vienna_to_ct(seq, pairs))


if __name__ == "__main__":
  main()
