import pytest

from rnapy.model.parse.sequence import (
    db_to_secondary,
    default_brackets,
    secondary_to_db,
    seq_to_primary,
)


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


def test_extended_dot_bracket_uses_uppercase_openers() -> None:
    secondary = [5, 6, 7, 8, 9, 0, 1, 2, 3, 4]

    db = secondary_to_db(secondary)

    assert db == "([{<A)]}>a"
    assert db_to_secondary(db) == secondary


@pytest.mark.parametrize("secondary", [[5, 7, 4, 6, 2, 0, 3, 1], [6, 3, 5, 1, 7, 2, 0, 4]])
def test_secondary_to_db_recolors_representable_pseudoknot(secondary: list[int]) -> None:
    db = secondary_to_db(secondary)

    assert set(db) <= set("".join(default_brackets()[:2]))
    assert db_to_secondary(db) == secondary


def test_secondary_to_db_rejects_structure_that_needs_more_bracket_types() -> None:
    pair_count = len(default_brackets()) + 1
    secondary = list(range(pair_count, 2 * pair_count)) + list(range(pair_count))

    with pytest.raises(ValueError, match="not enough bracket types"):
        secondary_to_db(secondary)


def test_secondary_to_db_rejects_self_pairing() -> None:
    with pytest.raises(ValueError):
        secondary_to_db([0])


def test_secondary_to_db_rejects_unpaired_closing() -> None:
    with pytest.raises(ValueError):
        secondary_to_db([-1, 0])


@pytest.mark.parametrize("secondary", [[2, -1], [1, -1], [-2]])
def test_secondary_to_db_rejects_invalid_pairs(secondary: list[int]) -> None:
    with pytest.raises(ValueError):
        secondary_to_db(secondary)
