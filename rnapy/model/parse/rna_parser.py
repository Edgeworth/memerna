# Copyright 2022 E.


import re
from collections import deque
from decimal import Decimal

from rnapy.model.parse.sequence import db_to_secondary, seq_to_primary
from rnapy.model.rna import Rna


class RnaParser:
    @staticmethod
    def _extract_energy(data: str) -> Decimal | None:
        pattern = r"ENERGY\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+))"
        match = re.search(pattern, data, re.IGNORECASE)
        if not match:
            return None
        return Decimal(match.group(1))

    @staticmethod
    def from_ct_file(data: str) -> Rna:
        lines = data.strip().splitlines()

        name = lines[0].strip().split(maxsplit=1)[-1]
        energy = RnaParser._extract_energy(name)
        line_fields = [re.split(r"\s+", i.strip()) for i in lines[1:] if i]
        seq = "".join(i[1] for i in line_fields)
        pairs = [-1 for i in range(len(seq))]
        for i, v in enumerate(line_fields):
            base = v[1]
            base_idx, prev_idx, next_idx, pair_idx = (int(v[0]), int(v[2]), int(v[3]), int(v[4]))
            if base_idx != i + 1:
                raise ValueError(f"base idx {base_idx} != {i + 1}")
            if base not in "GUAC":  # Only consider fully determined sequences for now.
                raise ValueError(f"base {base} not in GUAC")
            if prev_idx != base_idx - 1:
                raise ValueError(f"prev_idx {prev_idx} != {base_idx - 1}")
            if i < len(line_fields) - 1 and next_idx != base_idx + 1:
                raise ValueError(f"next_idx {next_idx} != {base_idx + 1}")
            if pair_idx != 0:
                pairs[pair_idx - 1] = base_idx - 1
                pairs[base_idx - 1] = pair_idx - 1
        return Rna(name=name, r=seq, s=pairs, energy=energy)

    @staticmethod
    def parse(
        *,
        name: str | None = None,
        seq: str | None = None,
        db: str | None = None,
        energy: str | None = None,
    ) -> Rna:
        return Rna(
            name=name,
            r=seq_to_primary(seq) if seq else None,
            s=db_to_secondary(db) if db else None,
            energy=Decimal(energy) if energy else None,
        )

    @staticmethod
    def from_db_file(data: str) -> Rna:
        name, seq, db = data.strip().splitlines()
        name, seq, db = name.strip(), seq.strip(), db.strip()
        name = re.sub(r"^> ", "", name)
        return RnaParser.parse(name=name, seq=seq, db=db)

    @staticmethod
    def from_any_file(data: str) -> Rna:
        if data[0] == ">":
            return RnaParser.from_db_file(data)
        return RnaParser.from_ct_file(data)

    @staticmethod
    def multi_from_ct_file(data: str) -> list[Rna]:
        q = deque(data.strip().splitlines())
        rnas = []
        while len(q) > 0:
            match = re.search(r"(\d+)", q[0].strip())
            if match is None:
                raise ValueError(f"Invalid CT file: {q[0]}")
            length = int(match.group(0))
            subdata = f"{q[0]}\n"
            q.popleft()
            for _ in range(length):
                subdata += f"{q[0]}\n"
                q.popleft()
            rnas.append(RnaParser.from_ct_file(subdata))
        return rnas

    @staticmethod
    def to_ct_file(rna: Rna) -> str:
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")
        if rna.s is None:
            raise ValueError(f"RNA {rna.name} has no secondary structure")

        name = rna.name if rna.name else "unnamed"
        ct = [f"{len(rna.r)}\t{name}"]

        for i, v in enumerate(rna.r):
            ct.append(f"{i + 1}\t{v}\t{i}\t{i + 2}\t{rna.s[i] + 1}\t{i + 1}")

        return "\n".join(ct)

    @staticmethod
    def to_db_file(rna: Rna) -> str:
        return f"> {rna.name}\n{rna.r}\n{rna.db()}\n"

    # See http://rna.urmc.rochester.edu/Text/File_Formats.html for this format.
    @staticmethod
    def to_seq_file(rna: Rna) -> str:
        return f";\n{rna.name}\n{rna.r}1"
