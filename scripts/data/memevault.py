# Copyright 2016 Eliot Courtney.
import glob
from pathlib import Path
import random
import sqlite3

from scripts.model.rna import Rna


# TODO: review, fix
class MemeVault:
    def __init__(self, path: Path, dataset: str):
        self.db = sqlite3.connect(path)
        self.dataset = dataset

    def add(self, rna):
        self.db.execute(
            f"INSERT INTO {self.dataset} VALUES (?, ?, ?)",
            (rna.name, rna.r, rna.db()),
        )
        self.db.commit()

    def get_with_seq(self, seq):
        c = self.db.execute(f"SELECT *  FROM {self.dataset} WHERE seq=?", (seq,))
        name, seq, db = c.fetchone()
        return Rna.from_name_seq_db(name, seq, db)

    def add_in_dir(self, dir):
        dir = fix_path(dir)
        for filename in glob.iglob(os.path.join(dir, "*.ct"), recursive=True):
            name = os.path.splitext(os.path.basename(filename))[0]
            rna = Rna.from_any_file(read_file(filename))
            rna.name = name
            if rna in self:
                print(f"Found duplicate RNAs: {rna.name}")
            self.add(rna)

    def add_random(self, l):
        seq = "".join(random.choice("GUAC") for _ in range(l))
        rna = Rna(name=f"len_{l}", r=seq, s=[-1] * l)
        self.add(rna)

    def __contains__(self, item):
        c = self.db.execute(f"SELECT count(*) FROM {self.dataset} WHERE name=?", (item.name,))
        return c.fetchone()[0] == 1

    def __getitem__(self, item):
        c = self.db.execute(f"SELECT * FROM {self.dataset} WHERE name=?", (item,))
        name, seq, db = c.fetchone()
        return Rna.from_name_seq_db(name, seq, db)

    def __iter__(self):
        c = self.db.execute(f"SELECT * FROM {self.dataset}")
        for name, seq, db in c.fetchall():
            yield Rna.from_name_seq_db(name, seq, db)
