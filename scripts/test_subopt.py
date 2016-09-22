#!/usr/bin/env python3
from common import *

icam1 = 'GCGCCCCAGUCGACGCUGAGCUCCUCUGCUACUCAGAGUUGCAACCUCAGCCUCGCUAUGGCUCCCAGCAGCCCCCGGCCCGCGCUGCCCGCACUCCUGGUCCUGCUCGGGGCUCUGUUCCCAGGACCUGGCAAUGCCCAGACAUCUGUGUCCCCCUCAAAAGUCAUCCUGCCCCGGGGAGGCUCCGUGCUGGUGACAUGCAGCACCUCCUGUGACCAGCCCAAGUUGUUGGGCAUAGAGACCCCGUUGCCUAAAAAGGAGUUGCUCCUGCCUGGGAACAACCGGAAGGUGUAUGAACUGAGCAAUGUGCAAGAAGAUAGCCAACCAAUGUGCUAUUCAAACUGCCCUGAUGGGCAGUCAACAGCUAAAACCUUCCUCACCGUGUACUGGACUCCAGAACGGGUGGAACUGGCACCCCUCCCCUCUUGGCAGCCAGUGGGCAAGAACCUUACCCUACGCUGCCAGGUGGAGGGUGGGGCACCCCGGGCCAACCUCACCGUGGUGCUGCUCCGUGGGGAGAAGGAGCUGAAACGGGAGCCAGCUGUGGGGGAGCCCGCUGAGGUCACGACCACGGUGCUGGUGAGGAGAGAUCACCAUGGAGCCAAUUUCUCGUGCCGCACUGAACUGGACCUGCGGCCCCAAGGGCUGGAGCUGUUUGAGAACACCUCGGCCCCCUACCAGCUCCAGACCUUUGUCCUGCCAGCGACUCCCCCACAACUUGUCAGCCCCCGGGUCCUAGAGGUGGACACGCAGGGGACCGUGGUCUGUUCCCUGGACGGGCUGUUCCCAGUCUCGGAGGCCCAGGUCCACCUGGCACUGGGGGACCAGAGGUUGAACCCCACAGUCACCUAUGGCAACGACUCCUUCUCGGCCAAGGCCUCAGUCAGUGUGACCGCAGAGGACGAGGGCACCCAGCGGCUGACGUGUGCAGUAAUACUGGGGAACCAGAGCCAGGAGACACUGCAGACAGUGACCAUCUACAGCUUUCCGGCGCCCAACGUGAUUCUGACGAAGCCAGAGGUCUCAGAAGGGACCGAGGUGACAGUGAAGUGUGAGGCCCACCCUAGAGCCAAGGUGACGCUGAAUGGGGUUCCAGCCCAGCCACUGGGCCCGAGGGCCCAGCUCCUGCUGAAGGCCACCCCAGAGGACAACGGGCGCAGCUUCUCCUGCUCUGCAACCCUGGAGGUGGCCGGCCAGCUUAUACACAAGAACCAGACCCGGGAGCUUCGUGUCCUGUAUGGCCCCCGACUGGACGAGAGGGAUUGUCCGGGAAACUGGACGUGGCCAGAAAAUUCCCAGCAGACUCCAAUGUGCCAGGCUUGGGGGAACCCAUUGCCCGAGCUCAAGUGUCUAAAGGAUGGCACUUUCCCACUGCCCAUCGGGGAAUCAGUGACUGUCACUCGAGAUCUUGAGGGCACCUACCUCUGUCGGGCCAGGAGCACUCAAGGGGAGGUCACCCGCGAGGUGACCGUGAAUGUGCUCUCCCCCCGGUAUGAGAUUGUCAUCAUCACUGUGGUAGCAGCCGCAGUCAUAAUGGGCACUGCAGGCCUCAGCACGUACCUCUAUAACCGCCAGCGGAAGAUCAAGAAAUACAGACUACAACAGGCCCAAAAAGGGACCCCCAUGAAACCGAACACACAAGCCACGCCUCCCUGAACCUAUCCCGGGACAGGGCCUCUUCCUCGGCCUUCCCAUAUUGGUGGCAGUGGUGCCACACUGAACAGAGUGGAAGACAUAUGCCAUGCAGCUACACCUACCGGCCCUGGGACGCCGGAGGACAGGGCAUUGUCCUCAGUCAGAUACAACAGCAUUUGGGGCCAUGGUACCUGCACACCUAAAACACUAGGCCACGCAUCUGAUCUGUAGUCACAUGACUAAGCCAAGAGGAAGGAGCAAGACUCAAGACAUGAUUGAUGGAUGUUAAAGUCUAGCCUGAUGAGAGGGGAAGUGGUGGGGGAGACAUAGCCCCACCAUGAGGACAUACAACUGGGAAAUACUGAAACUUGCUGCCUAUUGGGUAUGCUGAGGCCCACAGACUUACAGAAGAAGUGGCCCUCCAUAGACAUGUGUAGCAUCAAAACACAAAGGCCCACACUUCCUGACGGAUGCCAGCUUGGGCACUGCUGUCUACUGACCCCAACCCUUGAUGAUAUGUAUUUAUUCAUUUGUUAUUUUACCAGCUAUUUAUUGAGUGUCUUUUAUGUAGGCUAAAUGAACAUAGGUCUCUGGCCUCACGGAGCUCCCAGUCCAUGUCACAUUCAAGGUCACCAGGUACAGUUGUACAGGUUGUACACUGCAGGAGAGUGCCUGGCAAAAAGAUCAAAUGGGGCUGGGACUUCUCAUUGGCCAACCUGCCUUUCCCCAGAAGGAGUGAUUUUUCUAUCGGCACAAAAGCACUAUAUGGACUGGUAAUGGUUCACAGGUUCAGAGAUUACCCAGUGAGGCCUUAUUCCUCCCUUCCCCCCAAAACUGACACCUUUGUUAGCCACCUCCCCACCCACAUACAUUUCUGCCAGUGUUCACAAUGACACUCAGCGGUCAUGUCUGGACAUGAGUGCCCAGGGAAUAUGCCCAAGCUAUGCCUUGUCCUCUUGUCCUGUUUGCAUUUCACUGGGAGCUUGCACUAUUGCAGCUCCAGUUUCCUGCAGUGAUCAGGGUCCUGCAAGCAGUGGGGAAGGGGGCCAAGGUAUUGGAGGACUCCCUCCCAGCUUUGGAAGGGUCAUCCGCGUGUGUGUGUGUGUGUAUGUGUAGACAAGCUCUCGCUCUGUCACCCAGGCUGGAGUGCAGUGGUGCAAUCAUGGUUCACUGCAGUCUUGACCUUUUGGGCUCAAGUGAUCCUCCCACCUCAGCCUCCUGAGUAGCUGGGACCAUAGGCUCACAACACCACACCUGGCAAAUUUGAUUUUUUUUUUUUUUUUCAGAGACGGGGUCUCGCAACAUUGCCCAGACUUCCUUUGUGUUAGUUAAUAAAGCUUUCUCAACUGCC'
lendeltanum = [
  (400, 10, 5000), (500, 10, 5000), (800, 5, 5000), (900, 5, 5000),
  (1000, 4, 100000), (1200, 4, 100000), (1400, 3, 100000), (1600, 3, 100000)]


def run_num(alg, ldn):
  print('alg%s - nums' % alg)
  for l, _, num in ldn:
    res = run_command(
      './build/subopt', '-q', '-num', str(num), '-subopt-alg', alg, icam1[:l])
    print('  len %d, num %d: %s' % (l, num, res))


def run_delta(alg, ldn):
  print('alg%s - deltas' % alg)
  for l, d, _ in ldn:
    res = run_command(
      './build/subopt', '-q', '-delta', str(d), '-subopt-alg', alg, icam1[:l])
    print('  len %d, delta %d: %s' % (l, d, res))

def main():
  # run_delta('0', lendeltanum[:4])
  # run_delta('1', lendeltanum)
  # run_num('0', lendeltanum[:4])
  # run_num('1', lendeltanum)
  for delta in range(1, 15):
    os.system('./scripts/run.py -b random -s %d -rh' % delta)
    os.system('./scripts/run.py -b random -s %d -rd' % delta)
    os.system('./scripts/run.py -b random -s %d -vd2' % delta)
    os.system('./scripts/run.py -b random -s %d -vd3' % delta)
    os.system('./scripts/run.py -b random -s %d -k' % delta)



if __name__ == '__main__':
  main()
