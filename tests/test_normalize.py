"""Circular origins are biological coordinates, rather than string offsets."""

import random

import pytest

from stump_of_life.normalize import normalize_to_reference


def _dna(length: int = 320, seed: int = 741) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def _rotate(sequence: str, offset: int) -> str:
    return sequence[offset:] + sequence[:offset]


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def test_identical_reference_special_case() -> None:
    result = normalize_to_reference("acgt", "ACGT")
    assert (result.sequence, result.strand, result.rotation, result.edit_distance) == (
        "ACGT", "+", 0, 0,
    )


@pytest.mark.parametrize("reverse", [False, True])
def test_every_rotation_is_exactly_normalized(reverse: bool) -> None:
    reference = _dna(123)
    oriented = _reverse_complement(reference) if reverse else reference
    for offset in range(len(reference)):
        query = _rotate(oriented, offset)
        result = normalize_to_reference(query, reference, kmer_size=11)
        assert result.sequence == reference
        assert result.strand == ("-" if reverse else "+")
        assert result.edit_distance == 0
        transformed = _reverse_complement(query) if result.strand == "-" else query
        assert _rotate(transformed, result.rotation) == result.sequence


@pytest.mark.parametrize("position", [0, 1, 22, 150, 319])
@pytest.mark.parametrize("reverse", [False, True])
def test_substitution_survives_normalization(position: int, reverse: bool) -> None:
    reference = _dna()
    alternate = "A" if reference[position] != "A" else "C"
    mutated = reference[:position] + alternate + reference[position + 1:]
    query = _rotate(mutated, 97)
    if reverse:
        query = _reverse_complement(query)
    result = normalize_to_reference(query, reference)
    assert result.sequence == mutated
    assert result.edit_distance == 1


@pytest.mark.parametrize("reverse", [False, True])
def test_insertion_before_reference_origin_is_trailing(reverse: bool) -> None:
    reference = _dna()
    inserted = "TTAAC"
    query = _rotate(inserted + reference, 101)
    if reverse:
        query = _reverse_complement(query)
    result = normalize_to_reference(query, reference)
    assert result.sequence == reference + inserted
    assert result.edit_distance == len(inserted)


@pytest.mark.parametrize("reverse", [False, True])
def test_deletion_across_old_origin_uses_shared_new_anchor(reverse: bool) -> None:
    reference = _dna()
    deleted = reference[2:-3]
    query = _rotate(deleted, 101)
    if reverse:
        query = _reverse_complement(query)
    result = normalize_to_reference(query, _rotate(reference, 70))
    assert result.sequence == reference[70:-3] + reference[2:70]
    assert result.edit_distance == 5


def test_deleted_origin_with_equal_distance_alternatives_is_rejected() -> None:
    reference = _dna()
    query = _rotate(reference[2:-3], 101)
    with pytest.raises(ValueError, match="origin is ambiguous"):
        normalize_to_reference(query, reference)


@pytest.mark.parametrize("position", [1, 8, 170, 317])
def test_internal_insertion_and_deletion_do_not_shift_origin(position: int) -> None:
    reference = _dna()
    inserted = reference[:position] + "GG" + reference[position:]
    result = normalize_to_reference(_rotate(inserted, 91), reference)
    assert result.sequence == inserted
    assert result.edit_distance == 2
    deleted = reference[:position] + reference[position + 2:]
    result = normalize_to_reference(_reverse_complement(_rotate(deleted, 29)), reference)
    assert result.sequence == deleted
    assert result.edit_distance == 2


@pytest.mark.parametrize("position", [0, 100, 319])
def test_reference_n_does_not_create_a_substitution(position: int) -> None:
    query = _dna()
    reference = query[:position] + "N" + query[position + 1:]
    result = normalize_to_reference(_rotate(query, 73), reference)
    assert result.sequence == query
    assert result.edit_distance == 0


def test_unknown_query_bases_are_preserved() -> None:
    reference = _dna()
    query = reference[:45] + "N" * 20 + reference[65:]
    result = normalize_to_reference(_reverse_complement(_rotate(query, 73)), reference)
    assert result.sequence == query
    assert result.edit_distance == 0


def test_compatible_iupac_bases_are_preserved() -> None:
    reference = _dna()
    ambiguity = {"A": "R", "C": "Y", "G": "R", "T": "Y"}[reference[80]]
    query = reference[:80] + ambiguity + reference[81:]
    result = normalize_to_reference(_rotate(query, 43), reference)
    assert result.sequence == query
    assert result.edit_distance == 0


def test_repetitive_reference_requires_an_unambiguous_anchor() -> None:
    reference = "ACGT" * 80
    with pytest.raises(ValueError, match="repetitive|ambiguous"):
        normalize_to_reference(_rotate(reference, 1), reference)


def test_reverse_complement_symmetric_circle_is_rejected() -> None:
    half = _dna(160)
    reference = half + _reverse_complement(half)
    with pytest.raises(ValueError, match="strand is ambiguous"):
        normalize_to_reference(_rotate(reference, 33), reference)


def test_unrelated_sequences_are_rejected() -> None:
    with pytest.raises(ValueError, match="too dissimilar"):
        normalize_to_reference(_dna(seed=111), _dna(seed=222))


def test_excessive_length_difference_is_rejected() -> None:
    reference = _dna()
    with pytest.raises(ValueError, match="too dissimilar in length"):
        normalize_to_reference(reference + "A" * 100, reference)


@pytest.mark.parametrize("query", ["", "ACGT-ACGT", "ACGTZACGT"])
def test_invalid_input_is_rejected(query: str) -> None:
    with pytest.raises(ValueError):
        normalize_to_reference(query, _dna())
