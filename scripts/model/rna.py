# Copyright 2016 Eliot Courtney.
from dataclasses import dataclass
import re

from scripts.model.parse import db_to_secondary
from scripts.model.parse import secondary_to_db
from scripts.model.parse import seq_to_primary


@dataclass
class Rna:
    name: str | None = None
    r: str | None = None
    s: list[int] | None = None
    energy: int | None = None

    def __post_init__(self):
        if self.r and self.s:
            assert len(self.r) == len(self.s)

    def __str__(self):
        res = ""
        if self.name:
            res += f"{self.name}\n"
        if self.r:
            res += f"{self.r}\n"
        if self.s:
            res += f"{self.db()}\n"
        if self.energy:
            res += f"{self.energy}\n"
        return res

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
        name = self.name if self.name else "unnamed"
        ct = [f"{len(self.r)}\t{name}"]

        for i, v in enumerate(self.r):
            ct.append(f"{i + 1}\t{v}\t{i}\t{i + 2}\t{self.s[i] + 1}\t{i + 1}")

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
            assert base in "GUAC"  # Only consider fully determined sequences for now.
            assert prev_idx == base_idx - 1
            if i < len(data) - 1:
                assert next_idx == base_idx + 1
            if pair_idx != 0:
                pairs[pair_idx - 1] = base_idx - 1
                pairs[base_idx - 1] = pair_idx - 1
        return Rna(name=name, r=seq, s=pairs)

    @staticmethod
    def parse(
        *,
        name: str | None = None,
        seq: str | None = None,
        db: str | None = None,
        energy: int | None = None,
    ):
        return Rna(
            name=name,
            r=seq_to_primary(seq) if seq else None,
            s=db_to_secondary(db) if db else None,
            energy=energy,
        )

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
