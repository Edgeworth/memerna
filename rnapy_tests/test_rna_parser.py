from decimal import Decimal

import pytest

from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.rna import Rna


def test_from_db_file_parses_basic() -> None:
    rna = RnaParser.from_db_file("> test\nGUAC\n(())\n")
    assert rna.name == "test"
    assert rna.r == "GUAC"
    assert rna.db() == "(())"


def test_from_db_file_parses_header_without_space() -> None:
    rna = RnaParser.from_db_file(">test\nGUAC\n(())\n")
    assert rna.name == "test"
    assert rna.r == "GUAC"
    assert rna.db() == "(())"


def test_from_db_file_strips_trailing_header_fields() -> None:
    rna = RnaParser.from_db_file(">test some description\nGUAC\n(())\n")
    assert rna.name == "test"
    assert rna.r == "GUAC"
    assert rna.db() == "(())"


def test_to_db_file_roundtrip_basic() -> None:
    rna = RnaParser.parse(name="test", seq="GUAC", db="(())")
    data = RnaParser.to_db_file(rna)
    parsed = RnaParser.from_db_file(data)
    assert parsed.name == "test"
    assert parsed.r == "GUAC"
    assert parsed.db() == "(())"


def test_rna_normalizes_sequences() -> None:
    assert Rna(r="gtac").r == "GUAC"


def test_rna_rejects_empty_sequence() -> None:
    with pytest.raises(ValueError, match="sequence must be non-empty"):
        Rna(r="")


def test_rna_rejects_empty_secondary() -> None:
    with pytest.raises(ValueError, match="secondary structure must be non-empty"):
        Rna(s=[])


def test_from_ct_file_parses_basic() -> None:
    ct = (
        "# File example.ct\n"
        "# External source: example\n"
        "4 test\n1 G 0 2 4 1\n2 U 1 3 3 2\n3 A 2 4 2 3\n4 C 3 0 1 4"
    )
    rna = RnaParser.from_ct_file(ct)
    assert rna.name == "test"
    assert rna.r == "GUAC"
    assert rna.db() == "(())"


def test_from_ct_file_accepts_last_row_next_index_past_end() -> None:
    rna = RnaParser.from_ct_file("1 test\n1 A 0 2 0 1")

    assert rna.r == "A"
    assert rna.s == [-1]


def test_ct_file_roundtrips_energy() -> None:
    rna = Rna(name="test", r="A", s=[-1], energy=Decimal("-1.20"))

    parsed = RnaParser.from_ct_file(RnaParser.to_ct_file(rna))

    assert parsed.name == rna.name
    assert parsed.energy == rna.energy


def test_from_ct_file_rejects_wrong_record_length() -> None:
    ct = "4 test\n1 G 0 0 0 1"
    with pytest.raises(ValueError, match="declares 4 bases, found 1"):
        RnaParser.from_ct_file(ct)


def test_multi_from_ct_file_rejects_truncated_record() -> None:
    ct = "4 test\n1 G 0 0 0 1"
    with pytest.raises(ValueError, match="declares 4 bases, found 1"):
        RnaParser.multi_from_ct_file(ct)


def test_to_ct_file_last_row_next_is_zero() -> None:
    rna = Rna(name="ct", r="GUAC", s=[3, 2, 1, 0])
    ct = RnaParser.to_ct_file(rna)
    lines = ct.splitlines()
    assert rna.r is not None
    assert len(lines) == 1 + len(rna.r)
    fields = [line.split() for line in lines[1:]]
    assert len(fields[-1]) >= 4
    assert fields[-1][3] == "0"


def test_rna_rejects_malformed_secondary() -> None:
    with pytest.raises(ValueError, match="Invalid secondary structure"):
        Rna(name="invalid", r="GUAC", s=[3, -1, -1, -1])


def test_from_any_file_rejects_empty_input() -> None:
    with pytest.raises(ValueError):
        RnaParser.from_any_file("")


def test_from_any_file_ignores_leading_whitespace() -> None:
    rna = RnaParser.from_any_file("\n  > test\nGUAC\n(())\n")
    assert rna.name == "test"


def test_output_formats_require_their_rna_components() -> None:
    with pytest.raises(ValueError, match="name must be non-empty"):
        RnaParser.to_db_file(Rna(r="A", s=[-1]))
    with pytest.raises(ValueError, match="no secondary structure"):
        RnaParser.to_db_file(Rna(name="missing", r="A"))
    with pytest.raises(ValueError, match="no sequence"):
        RnaParser.to_seq_file(Rna(name="missing"))


def test_rna_rejects_mismatched_lengths() -> None:
    with pytest.raises(ValueError, match="same length"):
        Rna(r="GU", s=[-1])
