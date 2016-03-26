import re

MAX = 1000000000


def parse_stack_txt(data):
  """Outputs stacking data. Impossible is set to MAX. Output as integers ten times bigger than the original value.

  Outputs 16 matrices in AA, AC ... order. Each matrix is outputted in AA, AC ... order too.
  """
  lines = [i.strip() for i in re.sub(r' +', ' ', data).split('\n')]
  output = ''
  for i in range(len(lines)):
    if re.match(r'(\s*3\' <-- 5\'\s*){4}', lines[i]):
      matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 1:i + 5]]
      for m in range(4):
        for r in range(4):
          for c in range(4):
            val = matrix_lines[r][m * 4 + c]
            if val == '.':
              val = MAX
            else:
              val = int(val.replace('.', ''))
            output += '%d ' % val
          output += '\n'
        output += '\n'
      i += 5
  return output


stack_txt = parse_stack_txt(open('orig_data/stack.txt', encoding='utf-8').read())
print(stack_txt)
