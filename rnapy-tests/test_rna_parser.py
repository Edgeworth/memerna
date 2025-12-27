import pytest

from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.rna import Rna


def test_from_db_file_parses_basic() -> None:
    rna = RnaParser.from_db_file("> test\nGUAC\n(())\n")
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


def test_from_ct_file_parses_basic() -> None:
    ct = "4 test\n1 G 0 2 4 1\n2 U 1 3 3 2\n3 A 2 4 2 3\n4 C 3 0 1 4"
    rna = RnaParser.from_ct_file(ct)
    assert rna.name == "test"
    assert rna.r == "GUAC"
    assert rna.db() == "(())"


def test_to_ct_file_last_row_next_is_zero() -> None:
    rna = Rna(name="ct", r="GUAC", s=[3, 2, 1, 0])
    ct = RnaParser.to_ct_file(rna)
    lines = ct.splitlines()
    assert rna.r is not None
    assert len(lines) == 1 + len(rna.r)
    fields = [line.split() for line in lines[1:]]
    assert len(fields[-1]) >= 4
    assert fields[-1][3] == "0"


def test_from_any_file_rejects_empty_input() -> None:
    with pytest.raises(ValueError):
        RnaParser.from_any_file("")
