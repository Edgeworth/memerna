#!/usr/bin/env python3
import glob
import os
import re
import sqlite3
import urllib.request

from bs4 import BeautifulSoup
from rna import RNA
from common import *


class MemeVault:
  def __init__(self, dataset):
    self.db = sqlite3.connect('memevault/memevault.db')
    self.dataset = dataset

  def add(self, rna):
    self.db.execute(
      'INSERT INTO %s VALUES (?, ?, ?)' % self.dataset, (rna.name, rna.seq, rna.db()))
    self.db.commit()

  def get_with_seq(self, seq):
    c = self.db.execute('SELECT *  FROM %s WHERE seq=?' % self.dataset, (seq,))
    name, seq, db = c.fetchone()
    return RNA.from_name_seq_db(name, seq, db)

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

class RNAStrand:
  def __init__(self, data):
    self.rnas = []

    raw_data = re.findall(r'File (\S+).dp\n[^\n]+\n[^\n]+\n\n(.*?)\n\n', data, re.S)
    for name, data in raw_data:
      if name != 'SRP_00247': continue
      data = re.sub(r'\s', '', data)
      assert len(data) % 2 == 0
      self.rnas.append(
          RNA.from_name_seq_db(name, data[:len(data) // 2], data[len(data) // 2:]))

  def insert(self, memevault):
    for rna in self.rnas:
      if rna in memevault:
        print('Skipping %s' % rna)
        continue
      print('Scraping %s' % rna)
      self.scrape_rna(rna)
      memevault.add(rna)

  def trim_category(self, cat):
    return re.sub(r'[\[\]?:]', '', cat).strip()

  def scrape_rna(self, rna):
    response = urllib.request.urlopen(
      'http://www.rnasoft.ca/strand/show_results.php?molecule_ID=%s' % rna.name)
    soup = BeautifulSoup(response.read(), 'html.parser')
    table = soup.html.body.table.tr.td.table.tr.td.table.next_sibling
    rows = [i for i in table.find_all('tr') if i.td.next_sibling is not None]
    d = {self.trim_category(i.td.get_text()): i.td.next_sibling.get_text() for i in rows}
    assert rna.name == d['Molecule ID']

class ArchiveII:
  def __init__(self, dir):
    self.dir = fix_path(dir)

  def insert(self, memevault):
    print(self.dir)
    for filename in glob.iglob(os.path.join(self.dir, '*.ct'), recursive=True):
      name = os.path.splitext(os.path.basename(filename))[0]
      rna = RNA.from_any_file(read_file(filename))
      rna.name = name
      if rna in memevault:
         print('Skipping', rna)
      else:
        memevault.add(rna)

# memevault = MemeVault('archiveii')
# rnastrand_scraper = RNAStrand(open('scripts/rnastrand.data').read())
# rnastrand_scraper.insert(memevault)
#
# archiveii = ArchiveII('~/software/rna/rnark/RNASets/RNAStructureArchive')
# archiveii.insert(memevault)
