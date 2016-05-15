import re
import string


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

class RNA:
  def __init__(self, name=None, seq=None, pairs=None):
    self.name = name
    self.seq = seq
    self.pairs = pairs

  def __str__(self, *args, **kwargs):
    return '%s' % self.name

  # Handles pseudoknots.
  def db(self):
    brackets = ['()', '[]', '{}', '<>']
    brackets += ['%s%s' % i for i in zip(string.ascii_lowercase, string.ascii_uppercase)]
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
        db += brackets[stack_map[sid]][bid]
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
    # Only consider fully determined sequences with no pseudoknots for now.
    assert all(i in 'GUAC' for i in seq)
    assert all(i in '().' for i in db)
    stack = []
    pairs = [-1 for i in range(len(seq))]
    for i, v in enumerate(db):
      if v == '(':
        stack.append(i)
      elif v == ')':
        pair = stack.pop()
        pairs[pair] = i
        pairs[i] = pair
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
