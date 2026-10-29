"""The shipped real sequences must remain unaltered and origin-independent."""

import hashlib
import json
import os
import shutil
from pathlib import Path

import pytest

from stump_of_life.alignment import align_circular, reverse_complement
from stump_of_life.io import read_fasta

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "taiwan"


def test_public_sequence_checksums_and_single_ancient_reference() -> None:
    provenance = json.loads((EXAMPLE / "provenance.json").read_text())
    records = read_fasta(EXAMPLE / "mitogenomes.fasta")
    assert len(records) == 10
    assert len(records["KF540505.1"]) == 16569
    assert "N" not in records["KF540505.1"]
    for record in provenance["sequences"]:
        assert hashlib.sha256(records[record["accession"]].encode()).hexdigest() == record["sequence_sha256"]


def test_real_mitogenomes_keep_variant_calls_after_rotation_and_reverse_complement(tmp_path: Path) -> None:
    mafft = shutil.which(os.environ.get("MAFFT_PATH", "mafft"))
    if mafft is None:
        pytest.skip("Install MAFFT for the real mitochondrial-data integration test")
    records = read_fasta(EXAMPLE / "mitogenomes.fasta")
    kwargs = {"mafft": mafft, "threads": 1, "mafft_binaries": os.environ.get("MAFFT_BINARIES")}
    original = align_circular(records, "KF540505.1", tmp_path / "original", **kwargs)
    transformed: dict[str, str] = {}
    for index, (identifier, sequence) in enumerate(records.items()):
        if identifier == "KF540505.1":
            transformed[identifier] = sequence
            continue
        cut = 7901 + index * 79
        rotated = sequence[cut:] + sequence[:cut]
        transformed[identifier] = reverse_complement(rotated) if index % 2 else rotated
    changed = align_circular(transformed, "KF540505.1", tmp_path / "rotated", **kwargs)
    assert original["counts_by_sample"] == changed["counts_by_sample"]
    assert read_fasta(tmp_path / "original" / "aligned.fasta", aligned=True) == read_fasta(
        tmp_path / "rotated" / "aligned.fasta", aligned=True,
    )
    for identifier in ("KF540524.1", "KF540588.1"):
        assert changed["counts_by_sample"][identifier]["substitution"] == 4
