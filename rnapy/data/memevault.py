# Copyright 2016 E.
import random
import sqlite3
from collections.abc import Generator
from pathlib import Path

from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.rna import Rna


class MemeVault:
    dataset: str

    def __init__(self, path: Path, dataset: str) -> None:
        self.db = sqlite3.connect(path)
        self.dataset = dataset

    def maybe_create(self) -> None:
        self.db.execute(
            f"CREATE TABLE IF NOT EXISTS {self.dataset} "
            "(name TEXT PRIMARY KEY, seq TEXT NOT NULL, db TEXT NOT NULL)"
        )
        self.db.commit()

    def clear(self) -> None:
        self.db.execute(f"DELETE FROM {self.dataset}")
        self.db.commit()

    def add(self, rna: Rna) -> None:
        self.db.execute(f"INSERT INTO {self.dataset} VALUES (?, ?, ?)", (rna.name, rna.r, rna.db()))
        self.db.commit()

    def get_with_seq(self, seq: str) -> Rna:
        c = self.db.execute(f"SELECT *  FROM {self.dataset} WHERE seq=?", (seq,))
        name, seq, db = c.fetchone()
        return RnaParser.parse(name=name, seq=seq, db=db)

    def add_in_dir(self, dir_path: Path) -> None:
        for path in dir_path.glob("**/*.ct"):
            rna = RnaParser.from_any_file(path.read_text())
            rna.name = path.stem
            if rna in self:
                print(f"Found duplicate RNAs: {rna.name}")
            self.add(rna)

    def add_random(self, length: int, name: str) -> None:
        seq = "".join(random.choice("GUAC") for _ in range(length))
        rna = Rna(name=name, r=seq, s=[-1] * length)
        self.add(rna)

    def __contains__(self, item: Rna) -> bool:
        c = self.db.execute(f"SELECT 1 FROM {self.dataset} WHERE name=?", (item.name,)).fetchone()
        return c is not None

    def __getitem__(self, item: str) -> Rna:
        c = self.db.execute(f"SELECT * FROM {self.dataset} WHERE name=?", (item,))
        name, seq, db = c.fetchone()
        return RnaParser.parse(name=name, seq=seq, db=db)

    def __iter__(self) -> Generator[Rna, None, None]:
        c = self.db.execute(f"SELECT * FROM {self.dataset}")
        for name, seq, db in c.fetchall():
            yield RnaParser.parse(name=name, seq=seq, db=db)
