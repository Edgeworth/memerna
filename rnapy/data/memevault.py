# Copyright 2016 Eliot Courtney.
import random
import sqlite3
from collections.abc import Generator
from pathlib import Path

from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.parse.sequence import seq_to_primary
from rnapy.model.rna import Rna


class MemeVault:
    dataset: str

    def __init__(self, path: Path, dataset: str) -> None:
        if not dataset or "\x00" in dataset:
            raise ValueError("MemeVault dataset name must be non-empty and contain no NUL bytes")
        self.db = sqlite3.connect(path)
        self.dataset = dataset
        escaped_dataset = dataset.replace('"', '""')
        self._table = f'"{escaped_dataset}"'

    def dataset_exists(self) -> bool:
        row = self.db.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (self.dataset,)
        ).fetchone()
        return row is not None

    def maybe_create(self) -> None:
        self.db.execute(
            f"CREATE TABLE IF NOT EXISTS {self._table} "
            "(name TEXT NOT NULL PRIMARY KEY, seq TEXT NOT NULL, db TEXT)"
        )
        self.db.commit()

    def clear(self) -> None:
        self.db.execute(f"DELETE FROM {self._table}")
        self.db.commit()

    def add(self, rna: Rna) -> None:
        name = rna.require_name()
        sequence = rna.require_sequence()
        db = rna.db() if rna.s is not None else None
        self.db.execute(f"INSERT INTO {self._table} VALUES (?, ?, ?)", (name, sequence, db))
        self.db.commit()

    def has_reference_structures(self) -> bool:
        row = self.db.execute(f"SELECT COUNT(*), COUNT(db) FROM {self._table}").fetchone()
        assert row is not None
        total, with_reference = int(row[0]), int(row[1])
        return total > 0 and with_reference == total

    def get_with_seq(self, seq: str) -> Rna:
        seq = seq_to_primary(seq)
        c = self.db.execute(f"SELECT * FROM {self._table} WHERE seq=?", (seq,))
        row = c.fetchone()
        if row is None:
            raise KeyError(seq)
        name, seq, db = row
        return RnaParser.parse(name=name, seq=seq, db=db)

    def add_in_dir(self, dir_path: Path) -> None:
        for path in dir_path.glob("**/*.ct"):
            rna = RnaParser.from_any_file(path.read_text())
            rna.name = path.stem
            if rna in self:
                print(f"Found duplicate RNAs: {rna.name}")
                continue
            self.add(rna)

    def add_random(self, length: int, name: str) -> None:
        seq = "".join(random.choice("GUAC") for _ in range(length))
        rna = Rna(name=name, r=seq)
        self.add(rna)

    def __contains__(self, item: Rna) -> bool:
        c = self.db.execute(f"SELECT 1 FROM {self._table} WHERE name=?", (item.name,)).fetchone()
        return c is not None

    def __getitem__(self, item: str) -> Rna:
        row = self.db.execute(f"SELECT * FROM {self._table} WHERE name=?", (item,)).fetchone()
        if row is None:
            raise KeyError(item)
        name, seq, db = row
        return RnaParser.parse(name=name, seq=seq, db=db)

    def __iter__(self) -> Generator[Rna]:
        c = self.db.execute(f"SELECT * FROM {self._table}")
        for name, seq, db in c.fetchall():
            yield RnaParser.parse(name=name, seq=seq, db=db)
