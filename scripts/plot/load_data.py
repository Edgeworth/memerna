import pandas as pd

from common import fix_path, read_file

colmap = {
  'name': 'Name',
  'length': 'Length (nuc)',
  'real': 'Wall time (s)',
  'usersys': 'User+sys time (s)',
  'maxrss': 'Maximum RSS (B)',
  'fscore': 'F-Score',
  'ppv': 'PPV',
  'sensitivity': 'Sensitivity',
  'mfe': 'MFE'
}


def read_fold_frame(filename, subset):
  cols = ['name', 'run', 'length', 'real', 'usersys', 'maxrss', 'fscore', 'ppv', 'sensitivity', 'mfe']
  frame = pd.read_csv(filename, delimiter=' ', header=None, names=cols)
  all_rows = []
  for name, group in frame.groupby('name'):
    assert len(group) == 5
    if name not in subset: continue
    # Assert that some columns are all the same
    for colname in ['length', 'fscore', 'ppv', 'sensitivity', 'mfe']:
      assert (group[colname] == group[colname].iloc[0]).all()
    rows = [[v for i, v in enumerate(group.iloc[0]) if i != 1] for _ in range(3)]
    # Remove outliers and average
    for colname in ['real', 'usersys', 'maxrss']:
      col = list(group[colname])
      col.sort()
      col = col[1:-1]
      for i in range(3):
        rows[i][cols.index(colname) - 1] = col[i]
    all_rows.extend(rows)
  cleaned_frame = pd.DataFrame(all_rows, columns=cols[:1] + cols[2:])
  assert len(all_rows) // 3 == len(subset)
  return cleaned_frame


def read_subopt_frame(filename):
  pass


class DataSet:
  def __init__(self, name, fmap):
    self.name = name
    self.fmap = fmap


def read_fold_dataset(name, filename_map, subset_file):
  subset = set(i.strip() for i in read_file(fix_path(subset_file)).strip().splitlines())
  frames = {name: read_fold_frame(filename, subset) for name, filename in filename_map.items()}
  return DataSet(name, frames)
