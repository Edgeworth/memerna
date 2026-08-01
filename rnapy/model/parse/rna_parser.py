# Copyright 2022 Eliot Courtney.


import re
from collections import deque
from decimal import Decimal, InvalidOperation

from rnapy.model.parse.sequence import db_to_secondary
from rnapy.model.rna import Rna

_CT_ENERGY_RE = re.compile(
    r"^ENERGY\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)(?=\s|$)", re.IGNORECASE
)


def _ct_lines(data: str) -> list[str]:
    lines = data.strip().splitlines()
    first_record = 0
    while first_record < len(lines) and lines[first_record].lstrip().startswith("#"):
        first_record += 1
    return lines[first_record:]


class RnaParser:
    @staticmethod
    def _parse_ct_header(header: str) -> tuple[int, str]:
        parts = header.strip().split(maxsplit=1)
        if not parts:
            raise ValueError("Empty CT header")
        try:
            length = int(parts[0])
        except ValueError as exc:
            raise ValueError(f"Invalid CT header: {header}") from exc
        if length <= 0:
            raise ValueError(f"Invalid CT record length: {length}")
        return length, parts[1] if len(parts) > 1 else ""

    @staticmethod
    def _parse_ct_metadata(data: str) -> tuple[str, Decimal | None]:
        match = _CT_ENERGY_RE.match(data)
        if not match:
            return data, None
        return data[match.end() :].strip(), Decimal(match.group(1))

    @staticmethod
    def from_ct_file(data: str) -> Rna:
        lines = _ct_lines(data)
        if not lines:
            raise ValueError("Empty CT file")

        length, metadata = RnaParser._parse_ct_header(lines[0])
        if len(lines) - 1 != length:
            raise ValueError(f"CT header declares {length} bases, found {len(lines) - 1}")
        name, energy = RnaParser._parse_ct_metadata(metadata)
        line_fields = [re.split(r"\s+", line.strip()) for line in lines[1:]]
        if any(len(fields) < 5 for fields in line_fields):
            raise ValueError("CT base rows must contain at least five fields")
        bases = []
        pairs = [-1] * length
        for i, v in enumerate(line_fields):
            base = v[1]
            try:
                base_idx, prev_idx, next_idx, pair_idx = (
                    int(v[0]),
                    int(v[2]),
                    int(v[3]),
                    int(v[4]),
                )
            except ValueError as exc:
                raise ValueError(f"Invalid integer field in CT base row {i + 1}") from exc
            if base_idx != i + 1:
                raise ValueError(f"base idx {base_idx} != {i + 1}")
            bases.append(base)
            if prev_idx != base_idx - 1:
                raise ValueError(f"prev_idx {prev_idx} != {base_idx - 1}")
            if i < len(line_fields) - 1 and next_idx != base_idx + 1:
                raise ValueError(f"next_idx {next_idx} != {base_idx + 1}")
            if i == len(line_fields) - 1 and next_idx not in (0, base_idx + 1):
                raise ValueError(f"next_idx {next_idx} must be 0 or {base_idx + 1}")
            pairs[i] = pair_idx - 1 if pair_idx else -1
        return Rna(name=name, r="".join(bases), s=pairs, energy=energy)

    @staticmethod
    def parse(
        *,
        name: str | None = None,
        seq: str | None = None,
        db: str | None = None,
        energy: str | None = None,
    ) -> Rna:
        try:
            parsed_energy = Decimal(energy) if energy is not None else None
        except InvalidOperation as exc:
            raise ValueError(f"Invalid RNA energy: {energy}") from exc
        return Rna(
            name=name,
            r=seq,
            s=db_to_secondary(db) if db is not None else None,
            energy=parsed_energy,
        )

    @staticmethod
    def from_db_file(data: str) -> Rna:
        lines = [line.strip() for line in data.strip().splitlines() if line.strip()]
        if len(lines) != 3:
            raise ValueError("DB file must have exactly 3 non-empty lines (header, seq, db)")

        header, seq, db = lines[0], lines[1], lines[2]
        if not header.startswith(">"):
            raise ValueError("DB file header must start with '>'")

        name = re.sub(r"^>\s*", "", header).strip()
        name = name.split(maxsplit=1)[0] if name else ""
        if not name:
            raise ValueError("DB file identifier must be non-empty")
        return RnaParser.parse(name=name, seq=seq, db=db)

    @staticmethod
    def from_any_file(data: str) -> Rna:
        data = data.lstrip()
        if not data:
            raise ValueError("Empty input")
        if data[0] == ">":
            return RnaParser.from_db_file(data)
        return RnaParser.from_ct_file(data)

    @staticmethod
    def multi_from_ct_file(data: str) -> list[Rna]:
        q = deque(_ct_lines(data))
        rnas = []
        while len(q) > 0:
            length, _ = RnaParser._parse_ct_header(q[0])
            if len(q) - 1 < length:
                raise ValueError(f"CT header declares {length} bases, found {len(q) - 1}")
            subdata = f"{q[0]}\n"
            q.popleft()
            for _ in range(length):
                subdata += f"{q[0]}\n"
                q.popleft()
            rnas.append(RnaParser.from_ct_file(subdata))
        return rnas

    @staticmethod
    def to_ct_file(rna: Rna) -> str:
        sequence = rna.require_sequence()
        secondary = rna.require_secondary()
        name = rna.name or "unnamed"
        if name.splitlines() != [name]:
            raise ValueError("CT RNA name must not contain line breaks")
        if name != name.strip():
            raise ValueError("CT RNA name must not have leading or trailing whitespace")
        if _CT_ENERGY_RE.match(name):
            raise ValueError("CT RNA name must not start with ENERGY metadata")
        if rna.energy is not None:
            name = f"ENERGY = {rna.energy}  {name}"
        ct = [f"{len(sequence)}\t{name}"]

        for i, base in enumerate(sequence):
            next_idx = i + 2 if i < len(sequence) - 1 else 0
            ct.append(f"{i + 1}\t{base}\t{i}\t{next_idx}\t{secondary[i] + 1}\t{i + 1}")

        return "\n".join(ct)

    @staticmethod
    def to_db_file(rna: Rna) -> str:
        sequence = rna.require_sequence()
        name = rna.require_name()
        if any(char.isspace() for char in name):
            raise ValueError("DB RNA name must not contain whitespace")
        return f"> {name}\n{sequence}\n{rna.db()}\n"

    # See http://rna.urmc.rochester.edu/Text/File_Formats.html for this format.
    @staticmethod
    def to_seq_file(rna: Rna) -> str:
        sequence = rna.require_sequence()
        name = rna.name or "unnamed"
        if name.splitlines() != [name]:
            raise ValueError("SEQ RNA name must not contain line breaks")
        return f";\n{name}\n{sequence}1\n"
