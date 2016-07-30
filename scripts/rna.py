import random
import re
import string

BRACKETS = ['()', '[]', '{}', '<>']
BRACKETS += ['%s%s' % i for i in zip(string.ascii_lowercase, string.ascii_uppercase)]


class RNAAccuracy:
  def __init__(self, ppv, sensitivity):
    self.ppv = ppv
    self.sensitivity = sensitivity
    self.fscore = 2.0 * sensitivity * ppv
    if self.fscore != 0:
      self.fscore /= ppv + sensitivity

  @staticmethod
  def from_rna(true, predicted):
    num_predicted = sum(1 for i in predicted.pairs if i != -1) / 2
    num_true = sum(1 for i in true.pairs if i != -1) / 2
    num_correct = 0
    for i, pair in enumerate(predicted.pairs):
      if pair != -1 and true.pairs[i] == pair:
        num_correct += 1
    num_correct /= 2
    assert num_correct <= num_predicted and num_correct <= num_true
    ppv, sensitivity = num_correct, num_correct
    if num_correct != 0:
      ppv /= num_predicted
      sensitivity /= num_true
    return RNAAccuracy(ppv=ppv, sensitivity=sensitivity)

  def __str__(self):
    return 'F-Score: %.2f - PPV: %.2f - Sensitivity: %.2f' % (self.fscore, self.ppv, self.sensitivity)


def rnas_from_db_list(data):
  lines = [i for i in data.split('\n') if i]
  rnas = []
  assert len(lines) % 3 == 0
  for i in range(len(lines) // 3):
    rnas.append(RNA.from_name_seq_db(lines[3 * i], lines[3 * i + 1], lines[3 * i + 2]))
  return rnas


def can_pair(a, b):
  if a == 'G':
    return b == 'C' or b == 'U'
  if a == 'U':
    return b == 'G' or b == 'A'
  if a == 'A':
    return b == 'U'
  if a == 'C':
    return b == 'G'


def generate_foldings(rna, pick_one):
  seq = rna.seq
  assert all(i in 'GUAC' for i in seq)
  rnas = []
  last = []
  pairs = [-1] * len(seq)
  idx = 0

  def gaf_internal(i):
    nonlocal last, idx
    if i == len(seq):
      if len(last) == 0:
        rnas.append(RNA(str(idx), seq, pairs[:]))
        idx += 1
        return True
      return False
    choices = [1, 2, 3]
    if pick_one:
      random.shuffle(choices)
    succeeded = False
    for choice in choices:
      if choice == 1:
        # 1. Try opening a new bracket.
        if len(seq) - i - 1 >= len(last) + 1:
          last.append(i)
          succeeded = gaf_internal(i + 1)
          last.pop(-1)
      elif choice == 2:
        # 2. Try closing a bracket.
        if last and i - last[-1] - 1 >= 3 and can_pair(seq[i], seq[last[-1]]):
          p = last.pop(-1)
          pairs[i] = p
          pairs[p] = i
          succeeded = gaf_internal(i + 1)
          pairs[i] = -1
          pairs[p] = -1
          last.append(p)
      else:
        # 3. Try doing nothing.
        if len(seq) - i - 1 >= len(last):
          succeeded = gaf_internal(i + 1)
      if succeeded and pick_one:
        break
    return succeeded

  gaf_internal(0)
  return rnas

# TODO: this is slow.
def generate_random_foldings(rna, num):
  rnas = []
  for i in range(num):
    rnas.extend(generate_foldings(rna, True))
  return rnas

def generate_all_foldings(rna):
  return generate_foldings(rna, False)


class RNA:
  def __init__(self, name=None, seq=None, pairs=None):
    self.name = name
    self.seq = seq
    self.pairs = pairs
    assert len(seq) == len(pairs)

  def __str__(self, *args, **kwargs):
    return '%s:\n  %s\n  %s' % (self.name, self.seq, self.db())

  # Handles pseudoknots.
  def db(self):
    inter = []
    popularity = {}
    stacks = []
    for i, pair in enumerate(self.pairs):
      assert pair != i
      if pair == -1:
        inter.append((-1, 0))
      elif pair > i:
        best = len(stacks)
        for sid, stack in enumerate(stacks):
          if not stack or pair < stack[-1]:
            best = sid
        if best == len(stacks):
          stacks.append([])
        stacks[best].append(pair)
        inter.append((best, 0))
      else:
        best = -1
        for sid, stack in enumerate(stacks):
          if stack and stack[-1] == i:
            best = sid
        assert best != -1
        stacks[best].pop()
        inter.append((best, 1))
        popularity.setdefault(best, 0)
        popularity[best] += 1
    db = ''
    stacks_by_pop = sorted(popularity.keys(), key=lambda x: popularity[x], reverse=True)
    stack_map = {v: i for i, v in enumerate(stacks_by_pop)}
    for sid, bid in inter:
      if sid == -1:
        db += '.'
      else:
        db += BRACKETS[stack_map[sid]][bid]
    return db

  def to_ct_file(self):
    ct = ['%d %s' % (len(self.seq), self.name)]

    for i, v in enumerate(self.seq):
      ct.append('%d %s %d %d %d %d' % (i + 1, v, i, i + 2, self.pairs[i] + 1, i + 1))

    return '\n'.join(ct)

  def to_db_file(self):
    return '> %s\n%s\n%s\n' % (self.name, self.seq, self.db())

  # See http://rna.urmc.rochester.edu/Text/File_Formats.html for this format.
  def to_seq_file(self):
    return ';\n%s\n%s1' % (self.name, self.seq)

  @staticmethod
  def from_ct_file(data):
    data = data.strip().split('\n')

    name = re.split(r'\s+', data[0].strip())[1]
    data = [re.split(r'\s+', i.strip()) for i in data[1:] if i]
    seq = ''.join(i[1] for i in data)
    pairs = [-1 for i in range(len(seq))]
    for i, v in enumerate(data):
      base = v[1]
      base_idx, prev_idx, next_idx, pair_idx = (
        int(v[0]), int(v[2]), int(v[3]), int(v[4]))
      assert base_idx == i + 1
      # Only consider fully determined sequences for now.
      assert base in "GUAC"
      assert prev_idx == base_idx - 1
      if i < len(data) - 1:
        assert next_idx == base_idx + 1
      if pair_idx != 0:
        pairs[pair_idx - 1] = base_idx - 1
        pairs[base_idx - 1] = pair_idx - 1
    return RNA(name=name, seq=seq, pairs=pairs)

  @staticmethod
  def from_name_seq_db(name, seq, db):
    # Only consider fully determined sequences.
    assert all(i in 'GUAC' for i in seq)
    opening = {v[0] for v in BRACKETS}
    closing = {v[1]: v[0] for v in BRACKETS}
    stack = {}
    pairs = [-1 for i in range(len(db))]
    for i, v in enumerate(db):
      if v in opening:
        stack.setdefault(v, [])
        stack[v].append(i)
      elif v in closing:
        pair = stack[closing[v]].pop()
        pairs[pair] = i
        pairs[i] = pair
      else:
        assert v == '.'
    return RNA(name=name, seq=seq, pairs=pairs)

  @staticmethod
  def from_db_file(data):
    name, seq, db = data.strip().split('\n')
    name, seq, db = name.strip(), seq.strip(), db.strip()
    name = re.sub(r'^> ', '', name)
    return RNA.from_name_seq_db(name, seq, db)

  @staticmethod
  def from_any_file(data):
    if data[0] == '>':
      return RNA.from_db_file(data)
    else:
      return RNA.from_ct_file(data)
