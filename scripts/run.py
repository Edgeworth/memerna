import argparse
import os
import shutil
import sys
import tempfile

from common import *
from rna import *
from scrape import MemeVault


class RNAstructure:
  def __init__(self, loc):
    self.loc = loc
    os.putenv('DATAPATH', os.path.join(loc, 'data_tables'))

  def fold(self, rna):
    with tempfile.NamedTemporaryFile('w') as f, tempfile.NamedTemporaryFile('r') as out:
      f.write(rna.to_seq_file())
      f.flush()
      benchmark_results, _ = benchmark_command(
        os.path.join(self.loc, 'exe', 'Fold'), '-mfe', f.name, out.name)
      predicted = RNA.from_any_file(out.read())
    return (predicted, benchmark_results)

  def efn(self, rna, logarithmic=False):
    with tempfile.NamedTemporaryFile('w') as f, tempfile.NamedTemporaryFile('r') as out:
      f.write(rna.to_ct_file())
      f.flush()
      extra_args = []
      if not logarithmic:
        extra_args.append('-s')
      run_command(os.path.join(self.loc, 'exe', 'efn2'), *extra_args, f.name, out.name)
      match = re.search('Energy = (.+)', out.read().strip())
      energy = float(match.group(1))
    return energy

  def close(self):
    pass

  def __str__(self):
    return 'RNAstructure'


# TODO: On hold until rnark works.
class Rnark:
  def __init__(self, loc):
    self.loc = loc

  def fold(self, rna):
    cwd = os.getcwd()
    os.chdir(os.path.dirname(rnark_loc))
    # run_command(os.path.join('.', os.path.basename(rnark_loc)), '1', rna.seq, rna.db)
    os.chdir(cwd)

  def efn(self, rna):
    cwd = os.getcwd()
    os.chdir(os.path.dirname(rnark_loc))
    # run_command(os.path.join('.', os.path.basename(rnark_loc)), '1', seq, pairs)
    os.chdir(cwd)

  def close(self):
    pass

  def __str__(self):
    return 'Rnark'


class ViennaRNA:
  def __init__(self, loc):
    self.loc = loc

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
    return (predicted, benchmark_results)

  def efn(self, rna):
    print("Running ViennaRNA")

  def close(self):
    pass

  def __str__(self):
    return 'ViennaRNA'


class UNAFold:
  def __init__(self, loc):
    self.loc = loc
    self.tempdir = tempfile.mkdtemp()
    os.putenv('UNAFOLDDAT', os.path.join(loc, 'data'))

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
    return (predicted, benchmark_results)

  def efn(self, rna):
    pass

  def close(self):
    shutil.rmtree(self.tempdir)

  def __str__(self):
    return 'UNAFold'


class MemeRNA:
  def __init__(self, loc):
    self.loc = loc

  def fold(self, rna):
    pass

  def efn(self, rna):
    run_command(self.loc, seq, pairs)

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
        if idx <= 2738:
          print('Skipping %d' % idx)
          idx += 1
          continue
        results = []
        print('Running %s on #%d %s' % (program, idx, rna.name))
        for i in range(2):
          predicted, result = program.fold(rna)
          results.append(result)
        combined = BenchmarkResults.combine_benchmarks(results)
        accuracy = RNAAccuracy.from_rna(rna, predicted)
        energy = rnastructure.efn(predicted)
        f.write('%d %.5f %.5f %.5f %.5f %.5f %.5f %.2f\n' % (
          len(rna.seq), combined.real, combined.usersys, combined.maxrss,
          accuracy.fscore, accuracy.ppv, accuracy.sensitivity, energy
        ))
        idx += 1
  rnastructure.close()


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
  parser.add_argument('-b', '--benchmark', action='store_true')
  args = parser.parse_args(sys.argv[1:] + list(extra_args))

  if bool(args.predict) + bool(args.energy) + bool(args.benchmark) != 1:
    parser.error('Exactly one of --predict, --energy, or --benchmark is required.')

  programs = []

  if args.rnastructure_loc:
    programs.append(RNAstructure(fix_path(args.rnastructure_loc)))
  if args.rnark_loc:
    programs.append(Rnark(fix_path(args.rnark_loc)))
  if args.viennarna_loc:
    programs.append(ViennaRNA(fix_path(args.viennarna_loc)))
  if args.unafold_loc:
    programs.append(UNAFold(fix_path(args.unafold_loc)))
  if args.memerna_loc:
    programs.append(MemeRNA(fix_path(args.memerna_loc)))

  if args.benchmark:
    process_benchmark(programs, args)
  else:
    if bool(args.file) == bool(args.memevault):
      parser.error('Exactly one of --file or --memevault is required.')

    if args.file:
      rna = RNA.from_any_file(read_file(args.file))
    elif args.memevault:
      memevault = MemeVault('archiveii')
      rna = memevault[args.memevault]
    for program in programs:
      if args.predict:
        print('Folding with %s:\n%s' % (program, program.fold(rna)[0]))
      if args.energy:
        print('Energy with %s: %f' % (program, program.efn(rna)))

  for program in programs:
    program.close()


if __name__ == '__main__':
  process_command()
