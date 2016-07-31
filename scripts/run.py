import argparse
import os
import re
import shutil
import sys
import tempfile

from common import *
from rna import *
from scrape import MemeVault


class RNAstructure:
  def __init__(self, loc=None):
    try:
      import default_paths
      loc = loc or default_paths.MILESRNASTRUCTURE_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)
    os.putenv('DATAPATH', os.path.join(self.loc, 'data_tables'))

  def fold(self, rna):
    with tempfile.NamedTemporaryFile('w') as f, tempfile.NamedTemporaryFile('r') as out:
      prev_dir = os.getcwd()
      os.chdir(self.loc)
      f.write(rna.to_seq_file())
      f.flush()
      benchmark_results, _ = benchmark_command(
        os.path.join('build', 'Fold'), '-mfe', f.name, out.name)
      output = out.read()
      predicted = RNA.from_any_file(output)
      os.chdir(prev_dir)
    return predicted, benchmark_results, output

  def batch_efn(self, rnas):
    prev_dir = os.getcwd()
    os.chdir(self.loc)
    input = ''.join(rna.to_db_file() for rna in rnas)
    benchmark_results, stdout = benchmark_command(
      os.path.join('build', 'memefn'), input=input)
    matches = re.findall(r'^Energy: (.+)$', stdout, re.M)
    energies = [float(i) for i in matches]
    os.chdir(prev_dir)
    return energies, benchmark_results, stdout

  def efn(self, rna, logarithmic=False):
    with tempfile.NamedTemporaryFile('w') as f, tempfile.NamedTemporaryFile('r') as out:
      f.write(rna.to_ct_file())
      f.flush()
      extra_args = ['-w']
      # Note that not giving this flag doesn't make it logarithmic.
      # RNAstructure 5.8 adds the logarithmic and asymmetry models together in this case.
      # RNAstructure also uses a coefficient of -6 for the number of branches, rather than
      # the fitted -9.
      if not logarithmic:
        extra_args.append('-s')
      benchmark_results, _ = benchmark_command(
        os.path.join(self.loc, 'build', 'efn2'), *extra_args, f.name, out.name)
      print(_)
      output = out.read()
      match = re.search(r'[eE]nergy = (.+)', output.strip())
      energy = float(match.group(1))
    return energy, benchmark_results, output

  def close(self):
    pass

  def __str__(self):
    return 'RNAstructure'


# TODO: On hold until rnark works.
class Rnark:
  def __init__(self, loc=None):
    try:
      import default_paths
      loc = loc or default_paths.RNARK_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)

  def fold(self, rna):
    cwd = os.getcwd()
    os.chdir(os.path.dirname(rnark_loc))
    # run_command(os.path.join('.', os.path.basename(rnark_loc)), '1', rna.seq, rna.db)
    os.chdir(cwd)

  def efn(self, rna, logarithmic=False):
    cwd = os.getcwd()
    os.chdir(os.path.join(self.loc, 'bin'))
    flag = '0'
    if logarithmic:
      flag = '1'
    benchmark_results, stdout = benchmark_command(
      fix_path('Memefn'), flag, rna.seq, rna.db())
    os.chdir(cwd)
    match = re.search(r'Total energy: (.+)', stdout)
    energy = float(match.group(1)) / 10.0
    return energy, benchmark_results, stdout

  def close(self):
    pass

  def __str__(self):
    return 'Rnark'


class ViennaRNA:
  def __init__(self, loc=None):
    try:
      import default_paths
      loc = loc or default_paths.VIENNARNA_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)

  def fold(self, rna):
    with tempfile.NamedTemporaryFile('w') as f:
      f.write(rna.seq)
      f.flush()
      benchmark_results, stdout = benchmark_command(
        os.path.join(self.loc, 'src', 'bin', 'RNAfold'),
        #  '-d3',
        '--noPS', '-i', f.name)
      seq, db = stdout.strip().split('\n')
      db = db.split(' ')[0]
      predicted = RNA.from_name_seq_db(rna.name, seq.strip(), db.strip())
    return predicted, benchmark_results, stdout

  def efn(self, rna):
    print('Not implemented yet')

  def close(self):
    pass

  def __str__(self):
    return 'ViennaRNA'


class UNAFold:
  def __init__(self, loc=None):
    try:
      import default_paths
      loc = loc or default_paths.UNAFOLD_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)
    self.tempdir = tempfile.mkdtemp()
    os.putenv('UNAFOLDDAT', os.path.join(self.loc, 'data'))

  def fold(self, rna):
    prev_dir = os.getcwd()
    os.chdir(self.tempdir)
    with open(os.path.join(self.tempdir, 'rna.seq'), 'w') as f:
      f.write(rna.to_db_file())
      f.flush()
      benchmark_results, _ = benchmark_command(
        os.path.join(self.loc, 'src', 'hybrid-ss-min'), f.name)
      predicted = RNA.from_any_file(read_file(os.path.splitext(f.name)[0] + '.ct'))
    os.chdir(prev_dir)
    return predicted, benchmark_results, None

  def efn(self, rna):
    pass

  def close(self):
    shutil.rmtree(self.tempdir)

  def __str__(self):
    return 'UNAFold'


class MemeRNA:
  def __init__(self, loc=None):
    try:
      import default_paths
      loc = loc or default_paths.MEMERNA_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)

  def fold(self, rna):
    pass

  def batch_efn(self, rnas):
    prev_dir = os.getcwd()
    os.chdir(self.loc)
    input = ''.join(rna.to_db_file() for rna in rnas)
    benchmark_results, stdout = benchmark_command(
      os.path.join('build', 'batch_energy'), input=input)
    matches = re.findall(r'^Energy: (.+)$', stdout, re.M)
    energies = [float(i) / 10.0 for i in matches]
    os.chdir(prev_dir)
    return energies, benchmark_results, stdout

  def efn(self, rna):
    prev_dir = os.getcwd()
    os.chdir(self.loc)
    benchmark_results, stdout = benchmark_command(
      os.path.join('build', 'compute_energy'), rna.seq, rna.db())
    match = re.search(r'^Energy: (.+)$', stdout, re.M)
    energy = float(match.group(1)) / 10.0
    os.chdir(prev_dir)
    return energy, benchmark_results, stdout

  def close(self):
    pass

  def __str__(self):
    return 'MemeRNA'


def process_benchmark(programs, args):
  memevault = MemeVault('random')
  rnastructure = RNAstructure(fix_path('~/software/rna/rnastructure'))
  for program in programs:
    with open('%s_%s.results' % (program, memevault.dataset), 'w') as f:
      idx = 1
      for rna in memevault:
        results = []
        print('Running %s on #%d %s' % (program, idx, rna.name))
        for i in range(2):
          predicted, result, _ = program.fold(rna)
          results.append(result)
        combined = BenchmarkResults.combine_benchmarks(results)
        accuracy = RNAAccuracy.from_rna(rna, predicted)
        energy, _, _ = rnastructure.efn(predicted)
        f.write('%d %.5f %.5f %.5f %.5f %.5f %.5f %.2f\n' % (
          len(rna.seq), combined.real, combined.usersys, combined.maxrss,
          accuracy.fscore, accuracy.ppv, accuracy.sensitivity, energy
        ))
        idx += 1
  rnastructure.close()


def process_command(*extra_args):
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--file', type=str)
  parser.add_argument('-k', '--memevault', type=str)
  parser.add_argument('--rnastructure-loc')
  parser.add_argument('--rnark-loc')
  parser.add_argument('--viennarna-loc')
  parser.add_argument('--unafold-loc')
  parser.add_argument('--memerna-loc')
  parser.add_argument('-pr', '--rnastructure', action='store_true')
  parser.add_argument('-pm', '--rnark', action='store_true')
  parser.add_argument('-pv', '--viennarna', action='store_true')
  parser.add_argument('-pu', '--unafold', action='store_true')
  parser.add_argument('-pk', '--memerna', action='store_true')
  parser.add_argument('-p', '--predict', action='store_true')
  parser.add_argument('-e', '--energy', action='store_true')
  parser.add_argument('-b', '--benchmark', action='store_true')
  parser.add_argument('cmd', type=str, nargs='*')
  args = parser.parse_args(sys.argv[1:] + list(*extra_args))

  if bool(args.predict) + bool(args.energy) + bool(args.benchmark) != 1:
    parser.error('Exactly one of --predict, --energy, or --benchmark is required.')

  programs = []

  if args.rnastructure_loc or args.rnastructure:
    programs.append(RNAstructure(args.rnastructure_loc))
  if args.rnark_loc or args.rnark:
    programs.append(Rnark(args.rnark_loc))
  if args.viennarna_loc or args.viennarna:
    programs.append(ViennaRNA(args.viennarna_loc))
  if args.unafold_loc or args.unafold:
    programs.append(UNAFold(args.unafold_loc))
  if args.memerna_loc or args.memerna:
    programs.append(MemeRNA(args.memerna_loc))

  if args.benchmark:
    process_benchmark(programs, args)
  else:
    if args.file:
      rna = RNA.from_any_file(read_file(args.file))
    if args.memevault:
      memevault = MemeVault('archiveii')
      rna = memevault[args.memevault]


    for program in programs:
      if args.predict:
        if args.cmd:
          if len(args.cmd) != 1:
            parser.error('Direct specification requires one argument for prediction.')
          seq = args.cmd[0]
          rna = RNA('user', seq, [-1] * len(seq))
        frna, benchmark_results, output = program.fold(rna)
        print('Folding %s with %s: %s\n  %s\n%s' % (
          rna.name, program, frna.db(), benchmark_results, output))
      if args.energy:
        if args.cmd:
          if len(args.cmd) != 2:
            parser.error('Direct specification requires two arguments for efn.')
          rna = RNA.from_name_seq_db('user', *args.cmd)
        energy, benchmark_results, output = program.efn(rna)
        print('Energy of %s with %s: %f\n  %s\n%s' % (
          rna.name, program, energy, benchmark_results, output))

  for program in programs:
    program.close()


if __name__ == '__main__':
  process_command()
