# Copyright 2016 E.
from pathlib import Path
import random
import sqlite3
from typing import Generator

from scripts.model.rna import Rna


# TODO: review, fix
class MemeVault:
    def __init__(self, path: Path, dataset: str):
        self.db = sqlite3.connect(path)
        self.dataset = dataset

    def add(self, rna: Rna) -> None:
        self.db.execute(
            f"INSERT INTO {self.dataset} VALUES (?, ?, ?)",
            (rna.name, rna.r, rna.db()),
        )
        self.db.commit()

    def get_with_seq(self, seq: str) -> Rna:
        c = self.db.execute(f"SELECT *  FROM {self.dataset} WHERE seq=?", (seq,))
        name, seq, db = c.fetchone()
        return Rna.parse(name=name, seq=seq, db=db)

    def add_in_dir(self, dir_path: Path) -> None:
        for path in dir_path.glob("**/*.ct"):
            rna = Rna.from_any_file(path.read_text())
            rna.name = path.stem
            if rna in self:
                print(f"Found duplicate RNAs: {rna.name}")
            self.add(rna)

    def add_random(self, length: int) -> None:
        seq = "".join(random.choice("GUAC") for _ in range(length))
        rna = Rna(name=f"len_{length}", r=seq, s=[-1] * length)
        self.add(rna)

    def __contains__(self, item: Rna) -> bool:
        c = self.db.execute(f"SELECT 1 FROM {self.dataset} WHERE name=?", (item.name,)).fetchone()
        return c is not None

    def __getitem__(self, item: str) -> Rna:
        c = self.db.execute(f"SELECT * FROM {self.dataset} WHERE name=?", (item,))
        name, seq, db = c.fetchone()
        return Rna.parse(name=name, seq=seq, db=db)

    def __iter__(self) -> Generator[Rna, None, None]:
        c = self.db.execute(f"SELECT * FROM {self.dataset}")
        for name, seq, db in c.fetchall():
            yield Rna.parse(name=name, seq=seq, db=db)
