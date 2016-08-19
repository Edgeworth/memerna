#!/usr/bin/env python3
import argparse
import shutil
import tempfile

from common import *
from memevault import MemeVault
from rna import *


class HarnessFolder:
  def __init__(self, loc, name, flag):
    try:
      import default_paths
      loc = loc or default_paths.MEMERNA_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)
    self.name = name
    self.flag = flag

  def fold(self, rna):
    prev_dir = os.getcwd()
    os.chdir(self.loc)
    benchmark_results, stdout = benchmark_command(
      os.path.join('build', 'harness'), '-f', self.flag, input=rna.seq)
    os.chdir(prev_dir)
    lines = stdout.strip().split('\n')
    assert len(lines) == 2
    return RNA.from_name_seq_db(rna.name, rna.seq, lines[1]), benchmark_results

  def batch_efn(self, rnas):
    prev_dir = os.getcwd()
    os.chdir(self.loc)
    input = '\n'.join('%s\n%s' % (rna.seq, rna.db()) for rna in rnas)
    benchmark_results, stdout = benchmark_command(
      os.path.join('build', 'harness'), '-e', self.flag, input=input)
    os.chdir(prev_dir)
    energies = [float(i) / 10.0 for i in stdout.strip().split('\n')]
    return energies, benchmark_results

  def efn(self, rna):
    energies, benchmark_results = self.batch_efn([rna])
    return energies[0], benchmark_results

  def close(self):
    pass

  def __str__(self):
    return self.name


class RNAstructure(HarnessFolder):
  def __init__(self, loc=None):
    super().__init__(loc, 'RNAstructure', '-r')


class Rnark(HarnessFolder):
  def __init__(self, loc=None):
    super().__init__(loc, 'Rnark', '-m')


class MemeRNA(HarnessFolder):
  def __init__(self, loc=None):
    super().__init__(loc, 'MemeRNA', '-k')


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
        #'-d3',
        '--noPS', '-i', f.name)
      seq, db = stdout.strip().split('\n')
      db = db.split(' ')[0]
      predicted = RNA.from_name_seq_db(rna.name, seq.strip(), db.strip())
    return predicted, benchmark_results

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
    return predicted, benchmark_results

  def efn(self, rna):
    pass

  def close(self):
    shutil.rmtree(self.tempdir)

  def __str__(self):
    return 'UNAFold'


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
  parser.add_argument('-p', '--path', type=str)
  parser.add_argument('-kv', '--memevault', type=str)
  parser.add_argument('--viennarna-loc')
  parser.add_argument('--unafold-loc')
  parser.add_argument('--memerna-loc')
  parser.add_argument('-r', '--rnastructure', action='store_true')
  parser.add_argument('-m', '--rnark', action='store_true')
  parser.add_argument('-v', '--viennarna', action='store_true')
  parser.add_argument('-u', '--unafold', action='store_true')
  parser.add_argument('-k', '--memerna', action='store_true')
  parser.add_argument('-f', '--fold', action='store_true')
  parser.add_argument('-e', '--energy', action='store_true')
  parser.add_argument('-b', '--benchmark', action='store_true')
  parser.add_argument('cmd', type=str, nargs='*')
  args = parser.parse_args(sys.argv[1:] + list(*extra_args))

  if bool(args.fold) + bool(args.energy) + bool(args.benchmark) != 1:
    parser.error('Exactly one of --fold, --energy, or --benchmark is required.')

  programs = []

  if args.rnastructure:
    programs.append(RNAstructure(args.memerna_loc))
  if args.rnark:
    programs.append(Rnark(args.memerna_loc))
  if args.viennarna:
    programs.append(ViennaRNA(args.viennarna_loc))
  if args.unafold:
    programs.append(UNAFold(args.unafold_loc))
  if args.memerna:
    programs.append(MemeRNA(args.memerna_loc))

  if args.benchmark:
    process_benchmark(programs, args)
  else:
    if args.path:
      rna = RNA.from_any_file(read_file(args.path))
    if args.memevault:
      memevault = MemeVault('archiveii')
      rna = memevault[args.memevault]

    for program in programs:
      if args.fold:
        if args.cmd:
          if len(args.cmd) != 1:
            parser.error('Direct specification requires one argument for prediction.')
          seq = args.cmd[0]
          rna = RNA('user', seq, [-1] * len(seq))
        frna, benchmark_results = program.fold(rna)
        print('Folding %s with %s: %s\n  %s\n' % (
          rna.name, program, frna.db(), benchmark_results))
      if args.energy:
        if args.cmd:
          if len(args.cmd) != 2:
            parser.error('Direct specification requires two arguments for efn.')
          rna = RNA.from_name_seq_db('user', *args.cmd)
        energy, benchmark_results = program.efn(rna)
        print('Energy of %s with %s: %f\n  %s\n' % (
          rna.name, program, energy, benchmark_results))

  for program in programs:
    program.close()


if __name__ == '__main__':
  process_command()
