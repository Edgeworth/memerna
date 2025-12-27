import pytest

from rnapy.model.parse.sequence import db_to_secondary, secondary_to_db, seq_to_primary


def test_seq_to_primary_uppercases() -> None:
    assert seq_to_primary("guac") == "GUAC"


def test_seq_to_primary_rejects_invalid_bases() -> None:
    with pytest.raises(ValueError):
        seq_to_primary("GUX")


def test_db_to_secondary_parses_simple_pairs() -> None:
    assert db_to_secondary("(())") == [3, 2, 1, 0]
    assert db_to_secondary("[]") == [1, 0]


def test_db_to_secondary_rejects_unmatched_closing_bracket() -> None:
    with pytest.raises(ValueError):
        db_to_secondary(")")


def test_db_to_secondary_rejects_unmatched_opening_bracket() -> None:
    with pytest.raises(ValueError):
        db_to_secondary("(")


def test_db_to_secondary_rejects_invalid_character() -> None:
    with pytest.raises(ValueError):
        db_to_secondary("!")


def test_secondary_to_db_roundtrip_simple() -> None:
    assert secondary_to_db([3, 2, 1, 0]) == "(())"


def test_secondary_to_db_rejects_self_pairing() -> None:
    with pytest.raises(ValueError):
        secondary_to_db([0])


def test_secondary_to_db_rejects_unpaired_closing() -> None:
    with pytest.raises(ValueError):
        secondary_to_db([-1, 0])
