from pathlib import Path

import pytest

from stump_of_life.analysis import call_state, column_map, informative_columns, validate_alignment
from stump_of_life.io import read_fasta, write_fasta


def test_fasta_preserves_ambiguity_and_order(tmp_path: Path) -> None:
    records = {"ancient": "ACGTNRYSWKMBDHV", "modern": "TGCA"}
    path = tmp_path / "input.fa"
    write_fasta(records, path)
    assert read_fasta(path) == records


@pytest.mark.parametrize("text", [">x\nACGT\n>x\nACGT", ">x\n", "ACGT", ">x\nAC-GT", ">x\nACTZ"])
def test_invalid_fasta_fails(text: str, tmp_path: Path) -> None:
    path = tmp_path / "bad.fa"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_fasta(path)


def test_seam_insertion_and_internal_gap_coordinates() -> None:
    mapping = column_map("AC--GT-")
    assert mapping[2] == {
        "alignment_column": 3, "reference_position": None,
        "insertion_after_reference_position": 2, "insertion_rank": 1,
    }
    assert mapping[3]["insertion_rank"] == 2
    assert mapping[-1]["insertion_after_reference_position"] == 4


def test_uncertainty_is_not_a_substitution_or_variable_site() -> None:
    assert call_state("A", "N") == "unknown"
    assert call_state("R", "A") == "unknown"
    assert call_state("A", "-") == "deletion"
    assert call_state("-", "G") == "insertion"
    assert informative_columns({"a": "ACG-", "b": "NC-T"}) == [2, 3]


def test_imported_leading_and_trailing_insertion_addresses_are_unique() -> None:
    mapping = column_map("--ACG-")
    addresses = [
        (row["insertion_after_reference_position"], row["insertion_rank"])
        for row in mapping if row["reference_position"] is None
    ]
    assert addresses == [(3, 2), (3, 3), (3, 1)]


def test_all_gap_reference_is_rejected() -> None:
    with pytest.raises(ValueError, match="only gaps"):
        validate_alignment({"ref": "---", "query": "ACG"}, "ref")
