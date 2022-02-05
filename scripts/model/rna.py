# Copyright 2016 Eliot Courtney.
import re
from collections import deque
from typing import Optional

from scripts.model.parse import db_to_secondary, secondary_to_db, seq_to_primary


class Rna:
    name: Optional[str]
    r: Optional[str]
    s: Optional[list[int]]

    def __init__(self, name: str = None, r: str = None, s: list[int] = None):
        self.name = name
        self.r = r
        self.s = s
        if r and s:
            assert len(r) == len(s)

    def __str__(self):
        return f"{self.name}:\n  {self.r}\n  {self.db()}"

    # Used by command line parsing, for example.
    def __len__(self):
        if self.r:
            return len(self.r)
        if self.s:
            return len(self.s)
        return 0

    def db(self):
        if not self.s:
            return ""
        return secondary_to_db(self.s)

    def to_ct_file(self):
        name = self.name
        if not name:
            name = "unnamed"
        ct = [f"{len(self.r)}\t{name}"]

        for i, v in enumerate(self.r):
            ct.append(
                f"{int(i + 1)}\t{v}\t{int(i)}\t{int(i + 2)}\t{int(self.s[i] + 1)}\t{int(i + 1)}",
            )

        return "\n".join(ct)

    def to_db_file(self):
        return f"> {self.name}\n{self.r}\n{self.db()}\n"

    # See http://rna.urmc.rochester.edu/Text/File_Formats.html for this format.
    def to_seq_file(self):
        return f";\n{self.name}\n{self.r}1"

    @staticmethod
    def from_ct_file(data):
        data = data.strip().split("\n")

        name = re.split(r"\s+", data[0].strip())[1]
        data = [re.split(r"\s+", i.strip()) for i in data[1:] if i]
        seq = "".join(i[1] for i in data)
        pairs = [-1 for i in range(len(seq))]
        for i, v in enumerate(data):
            base = v[1]
            base_idx, prev_idx, next_idx, pair_idx = (int(v[0]), int(v[2]), int(v[3]), int(v[4]))
            assert base_idx == i + 1
            # Only consider fully determined sequences for now.
            assert base in "GUAC"
            assert prev_idx == base_idx - 1
            if i < len(data) - 1:
                assert next_idx == base_idx + 1
            if pair_idx != 0:
                pairs[pair_idx - 1] = base_idx - 1
                pairs[base_idx - 1] = pair_idx - 1
        return Rna(name=name, r=seq, s=pairs)

    @staticmethod
    def from_name_seq_db(name, seq, db):
        return Rna(name=name, r=seq_to_primary(seq), s=db_to_secondary(db))

    @staticmethod
    def from_db_file(data):
        name, seq, db = data.strip().split("\n")
        name, seq, db = name.strip(), seq.strip(), db.strip()
        name = re.sub(r"^> ", "", name)
        return Rna.from_name_seq_db(name, seq, db)

    @staticmethod
    def from_any_file(data):
        if data[0] == ">":
            return Rna.from_db_file(data)
        else:
            return Rna.from_ct_file(data)
