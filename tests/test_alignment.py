"""Exercise circular normalisation, a real MSA, and the exported coordinates.

Integration tests use MAFFT_PATH (or mafft on PATH) and optional MAFFT_BINARIES.
They skip when MAFFT is not installed; no machine-specific location is embedded.
"""

import csv
import json
import os
import random
import shutil
from collections.abc import Mapping
from pathlib import Path

import pytest

from stump_of_life.alignment import (
    align_circular,
    choose_common_anchor,
    restore_reference_origin,
)
from stump_of_life.io import read_fasta


def _dna(length: int = 640, seed: int = 19073) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def _rotate(sequence: str, offset: int) -> str:
    return sequence[offset:] + sequence[:offset]


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _assert_same_circle(actual: str, expected: str) -> None:
    """The final common-column cut may differ from the deposited query origin."""
    assert len(actual) == len(expected)
    assert actual in expected + expected


@pytest.fixture
def mafft_config() -> tuple[str, str | None]:
    configured = os.environ.get("MAFFT_PATH")
    executable = shutil.which(configured or "mafft")
    if executable is None:
        if configured:
            pytest.fail(f"Configured MAFFT_PATH is not executable: {configured}")
        pytest.skip("MAFFT is not installed; set MAFFT_PATH to run real alignment tests")
    return executable, os.environ.get("MAFFT_BINARIES")


def _run(
    records: Mapping[str, str], output: Path, mafft_config: tuple[str, str | None],
) -> dict[str, object]:
    executable, binaries = mafft_config
    return align_circular(
        records, "ref", output, mafft=executable, threads=1,
        mafft_binaries=binaries, kmer_size=23,
    )


def _check_exported_column_map(aligned: Mapping[str, str], output: Path) -> None:
    """Every output column must have a unique, reconstructable reference address."""
    rows = _csv_rows(output / "columns.csv")
    reference = aligned["ref"]
    assert len(rows) == len(reference)
    assert [int(row["alignment_column"]) for row in rows] == list(range(1, len(reference) + 1))
    observed_positions = []
    insertion_addresses = []
    for index, (base, row) in enumerate(zip(reference, rows, strict=True)):
        if base == "-":
            assert row["reference_position"] == ""
            insertion_addresses.append((
                int(row["insertion_after_reference_position"]), int(row["insertion_rank"]),
            ))
            before = reference[:index].replace("-", "")
            assert int(row["insertion_after_reference_position"]) == (
                len(before) or len(reference.replace("-", ""))
            )
        else:
            observed_positions.append(int(row["reference_position"]))
            assert row["insertion_after_reference_position"] == ""
            assert row["insertion_rank"] == ""
    assert observed_positions == list(range(1, len(reference.replace("-", "")) + 1))
    assert len(insertion_addresses) == len(set(insertion_addresses))


def test_restore_origin_rotates_whole_columns_including_insertions() -> None:
    # The reference was cut after AC. Gap columns must move together with all bases.
    aligned = {"ref": "G-TAC-", "query": "GGTACC"}
    restored = restore_reference_origin(aligned, "ref", 4, 2)
    assert restored == {"ref": "AC-G-T", "query": "ACCGGT"}
    original_columns = list(zip(*aligned.values(), strict=True))
    assert list(zip(*restored.values(), strict=True)) == original_columns[3:] + original_columns[:3]


@pytest.mark.parametrize("reference", ["ACGT" * 160, "N" * 640])
def test_repeated_or_unknown_only_anchor_is_rejected(reference: str) -> None:
    records = {"ref": reference, "query": _rotate(reference, 19)}
    with pytest.raises(ValueError, match="anchor|ambiguous|repetitive"):
        choose_common_anchor(records, "ref", 23)


def test_missing_reference_is_rejected_before_alignment() -> None:
    with pytest.raises(ValueError, match="Reference ID"):
        choose_common_anchor({"query": _dna()}, "ref", 23)


def test_real_alignment_is_invariant_to_rotation_and_strand(
    tmp_path: Path, mafft_config: tuple[str, str | None],
) -> None:
    reference = _dna()
    records = {"ref": reference}
    for offset in (0, 1, 217, len(reference) - 1):
        records[f"forward_{offset}"] = _rotate(reference, offset)
        records[f"reverse_{offset}"] = _reverse_complement(_rotate(reference, offset))
    output = tmp_path / "identical"
    summary = _run(records, output, mafft_config)
    aligned = read_fasta(output / "aligned.fasta", aligned=True)
    assert set(aligned) == set(records)
    assert all(sequence.replace("-", "") == reference for sequence in aligned.values())
    assert _csv_rows(output / "variants.csv") == []
    assert summary["samples"] == len(records)
    assert summary["reference_length"] == len(reference)
    for name in (
        "aligned.fasta", "normalized.fasta", "normalization.json", "columns.csv",
        "variants.csv", "summary.json", "mafft.log",
    ):
        assert (output / name).is_file()
    normalized = read_fasta(output / "normalized.fasta")
    for sequence in normalized.values():
        _assert_same_circle(sequence, reference)
    provenance = json.loads((output / "normalization.json").read_text())
    assert set(provenance["samples"]) == set(records)
    for identifier, transform in provenance["samples"].items():
        original = records[identifier]
        oriented = _reverse_complement(original) if transform["strand"] == "-" else original
        assert _rotate(oriented, transform["rotation"]) == normalized[identifier]
    _check_exported_column_map(aligned, output)


def test_real_alignment_calls_first_and_last_base_once_at_correct_coordinates(
    tmp_path: Path, mafft_config: tuple[str, str | None],
) -> None:
    reference = _dna()
    expected: dict[str, tuple[int, str, str]] = {}
    records = {"ref": reference}
    for index in (0, 1, len(reference) - 1):
        alternate = next(base for base in "ACGT" if base != reference[index])
        mutated = reference[:index] + alternate + reference[index + 1:]
        identifier = f"snp_{index + 1}"
        records[identifier] = _reverse_complement(_rotate(mutated, 317))
        expected[identifier] = (index + 1, reference[index], alternate)
    output = tmp_path / "seam_snps"
    _run(records, output, mafft_config)
    aligned = read_fasta(output / "aligned.fasta", aligned=True)
    rows = _csv_rows(output / "variants.csv")
    assert len(rows) == len(expected)
    assert {row["sample_id"] for row in rows} == set(expected)
    for row in rows:
        assert row["state"] == "substitution"
        assert (int(row["reference_position"]), row["reference_base"], row["sample_base"]) == (
            expected[row["sample_id"]]
        )
        column = int(row["alignment_column"]) - 1
        assert aligned["ref"][column] == row["reference_base"]
        assert aligned[row["sample_id"]][column] == row["sample_base"]
    _check_exported_column_map(aligned, output)


def test_real_alignment_retains_indels_and_maps_seam_insertions(
    tmp_path: Path, mafft_config: tuple[str, str | None],
) -> None:
    reference = _dna()
    expected = {
        "ref": reference,
        "insert": reference[:220] + "ACT" + reference[220:],
        "delete": reference[:391] + reference[395:],
        "seam_insert": reference + "GGTCA",
    }
    records = {
        "ref": reference,
        "insert": _rotate(expected["insert"], 317),
        "delete": _reverse_complement(_rotate(expected["delete"], 113)),
        "seam_insert": _reverse_complement(_rotate(expected["seam_insert"], 197)),
    }
    output = tmp_path / "indels"
    _run(records, output, mafft_config)
    aligned = read_fasta(output / "aligned.fasta", aligned=True)
    assert aligned["ref"].replace("-", "") == reference
    assert len({len(row) for row in aligned.values()}) == 1
    for identifier, sequence in aligned.items():
        _assert_same_circle(sequence.replace("-", ""), expected[identifier])
    rows = _csv_rows(output / "variants.csv")
    assert sum(row["sample_id"] == "insert" and row["state"] == "insertion" for row in rows) == 3
    assert sum(row["sample_id"] == "delete" and row["state"] == "deletion" for row in rows) == 4
    seam_insertions = [
        row for row in rows if row["sample_id"] == "seam_insert" and row["state"] == "insertion"
    ]
    assert len(seam_insertions) == 5
    assert all(int(row["insertion_after_reference_position"]) == len(reference) for row in seam_insertions)
    assert [int(row["insertion_rank"]) for row in seam_insertions] == list(range(1, 6))
    assert "".join(row["sample_base"] for row in seam_insertions) == "GGTCA"
    _check_exported_column_map(aligned, output)


def test_real_alignment_preserves_iupac_ambiguity_as_unknown_calls(
    tmp_path: Path, mafft_config: tuple[str, str | None],
) -> None:
    reference = _dna()
    replacements = {45: "N", 180: "R" if reference[180] in "AG" else "Y", 589: "N"}
    query = "".join(replacements.get(index, base) for index, base in enumerate(reference))
    output = tmp_path / "ambiguity"
    _run({"ref": reference, "ambiguous": _reverse_complement(_rotate(query, 79))}, output, mafft_config)
    aligned = read_fasta(output / "aligned.fasta", aligned=True)
    assert aligned["ambiguous"].replace("-", "") == query
    rows = _csv_rows(output / "variants.csv")
    assert len(rows) == len(replacements)
    assert all(row["state"] == "unknown" for row in rows)
    assert {int(row["reference_position"]): row["sample_base"] for row in rows} == {
        position + 1: base for position, base in replacements.items()
    }
    _check_exported_column_map(aligned, output)
