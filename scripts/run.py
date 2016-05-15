import argparse
import os
import tempfile

import sys
from rna import *
from common import *
from scrape import MemeVault


class RNAstructure:
  def __init__(self, loc):
    self.loc = loc
    os.putenv('DATAPATH', os.path.join(loc, 'data_tables'))

  def fold(self, rna):
    print("Running RNAstructure")
    with tempfile.NamedTemporaryFile('w') as f, tempfile.NamedTemporaryFile('r') as out:
      f.write(rna.to_seq_file())
      f.flush()
      run_command(os.path.join(self.loc, 'exe', 'Fold'), f.name, out.name, '-y', '-mfe')
      predicted = RNA.from_any_file(out.read())
      print('RNAstructure:\n%s\n%s\n%s' % (rna.db(), predicted.db(), RNAAccuracy.from_rna(rna, predicted)))

  def efn(self, rna):
    print("Running RNAstructure")
    with tempfile.NamedTemporaryFile() as f, tempfile.NamedTemporaryFile() as out:
      f.write(db_to_ct(seq, pairs).encode('UTF-8'))
      f.flush()
      run_command(os.path.join(self.loc, 'exe', 'efn2'), '-w', f.name, out.name)
      print('RNAstructure: %s' % out.read().strip().decode('UTF-8'))

class Rnark:
  def __init__(self, loc):
    self.loc = loc

  def fold(self, rna):
    cwd = os.getcwd()
    os.chdir(os.path.dirname(rnark_loc))
    #run_command(os.path.join('.', os.path.basename(rnark_loc)), '1', rna.seq, rna.db)
    os.chdir(cwd)

  def efn(self, rna):
    cwd = os.getcwd()
    os.chdir(os.path.dirname(rnark_loc))
    run_command(os.path.join('.', os.path.basename(rnark_loc)), '1', seq, pairs)
    os.chdir(cwd)

class ViennaRNA:
  def __init__(self, loc):
    self.loc = loc

  def fold(self, rna):
    print("Running ViennaRNA")
    with tempfile.NamedTemporaryFile('w') as f:
      f.write(rna.seq)
      f.flush()
      res = run_command(
        os.path.join(self.loc, 'src', 'bin', 'RNAfold'),
        '-d3', '--noPS', '-i', f.name)
      seq, db = res.stdout.decode('UTF-8').strip().split('\n')
      db = db.split(' ')[0]
      predicted = RNA.from_name_seq_db(rna.name, seq.strip(), db.strip())
      print('Vienna:\n%s\n%s\n%s' % (rna.db(), predicted.db(), RNAAccuracy.from_rna(rna, predicted)))

def efn(self, rna):
    print("Running ViennaRNA")

class UNAFold:
  def __init__(self, loc):
    self.loc = loc

  def fold(self, rna):
    pass

  def efn(self, rna):
    pass

class MemeRNA:
  def __init__(self, loc):
    self.loc = loc

  def fold(self, rna):
    pass

  def efn(self, rna):
    print("Running memerna")
    run_command(self.loc, seq, pairs)

def process_command(*extra_args):
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--file')
  parser.add_argument('-k', '--memevault')
  parser.add_argument('--rnastructure-loc')
  parser.add_argument('--rnark-loc')
  parser.add_argument('--viennarna-loc')
  parser.add_argument('--unafold-loc')
  parser.add_argument('--memerna-loc')
  parser.add_argument('-p', '--predict', action='store_true')
  parser.add_argument('-e', '--energy', action='store_true')
  args = parser.parse_args(sys.argv[1:] + list(extra_args))

  if args.file:
    rna = RNA.from_any_file(read_file(args.file))
  elif args.memevault:
    memevault = MemeVault('archiveii')
    rna = memevault[args.memevault]
  else:
    parser.error('One of --file or --memevault is required.')
  program = None

  if args.rnastructure_loc:
    program = RNAstructure(fix_path(args.rnastructure_loc))
  elif args.rnark_loc:
    program = Rnark(fix_path(args.rnark_loc))
  elif args.viennarna_loc:
    program = ViennaRNA(fix_path(args.viennarna_loc))
  elif args.unafold_loc:
    program = UNAFold(fix_path(args.unafold_loc))
  elif args.memerna_loc:
    program = MemeRNA(fix_path(args.memerna_loc))

  print(rna.seq, rna.db())
  if args.predict:
    print('Folding')
    program.fold(rna)
  if args.energy:
    program.efn(rna)

if __name__ == '__main__':
  process_command()
