# Copyright 2016, E.
#
# This file is part of memerna.
#
# memerna is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with memerna.
# If not, see <http://www.gnu.org/licenses/>.
import glob
import random
import sqlite3

from scripts.common import *
from scripts.rna import RNA


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
        print('Found duplicate RNAs: %s' % rna.name)
      self.add(rna)

  def add_random(self, l):
    seq = ''.join(random.choice('GUAC') for _ in range(l))
    rna = RNA(name='len_%d' % l, seq=seq, pairs=[-1] * l)
    self.add(rna)

  def __contains__(self, item):
    c = self.db.execute('SELECT count(*) FROM %s WHERE name=?' % self.dataset, (item.name,))
    return c.fetchone()[0] == 1

  def __getitem__(self, item):
    c = self.db.execute('SELECT * FROM %s WHERE name=?' % self.dataset, (item,))
    name, seq, db = c.fetchone()
    return RNA.from_name_seq_db(name, seq, db)

  def __iter__(self):
    c = self.db.execute('SELECT * FROM %s' % self.dataset)
    for name, seq, db in c.fetchall():
      yield RNA.from_name_seq_db(name, seq, db)
