#!/usr/bin/env python3
import glob
import random
import sqlite3

from common import *
from rna import RNA


class MemeVault:
  def __init__(self, dataset):
    self.db = sqlite3.connect('scripts/memevault.db')
    self.dataset = dataset

  def add(self, rna):
    self.db.execute(
      'INSERT INTO %s VALUES (?, ?, ?)' % self.dataset, (rna.name, rna.seq, rna.db()))
    self.db.commit()

  def get_with_seq(self, seq):
    c = self.db.execute('SELECT *  FROM %s WHERE seq=?' % self.dataset, (seq,))
    name, seq, db = c.fetchone()
    return RNA.from_name_seq_db(name, seq, db)

  def add_in_dir(self, dir):
    dir = fix_path(dir)
    for filename in glob.iglob(os.path.join(dir, '*.ct'), recursive=True):
      name = os.path.splitext(os.path.basename(filename))[0]
      rna = RNA.from_any_file(read_file(filename))
      rna.name = name
      if rna in self:
        print('Skipping', rna)
      else:
        self.add(rna)

  def add_random(self, num):
    for i in range(1, num + 1):
      seq = ''.join(random.choice('GUAC') for _ in range(i))
      rna = RNA(name='len_%d' % i, seq=seq, pairs=[-1] * i)
      self.add(rna)

  def __contains__(self, item):
    c = self.db.execute('SELECT count(*) FROM %s WHERE seq=?' % self.dataset, (item.seq,))
    return c.fetchone()[0] == 1

  def __getitem__(self, item):
    c = self.db.execute('SELECT * FROM %s WHERE name=?' % self.dataset, (item,))
    name, seq, db = c.fetchone()
    return RNA.from_name_seq_db(name, seq, db)

  def __iter__(self):
    c = self.db.execute('SELECT * FROM %s' % self.dataset)
    for name, seq, db in c.fetchall():
      yield RNA.from_name_seq_db(name, seq, db)
