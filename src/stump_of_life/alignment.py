"""Common circular anchor → conventional MSA → restore reference base one."""

import hashlib
import json
import os
import subprocess
from collections import Counter
from collections.abc import Mapping
from dataclasses import asdict
from pathlib import Path
from tempfile import TemporaryDirectory

from stump_of_life.analysis import validate_alignment, write_analysis
from stump_of_life.io import IUPAC, read_fasta, write_fasta
from stump_of_life.normalize import normalize_to_reference

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _kmers(sequence: str, size: int) -> list[str]:
    doubled = sequence + sequence[:size - 1]
    return [doubled[index:index + size] for index in range(len(sequence))]


def choose_common_anchor(
    records: Mapping[str, str], reference_id: str, kmer_size: int = 31,
) -> int:
    """Find an exact A/C/G/T anchor occurring once across both strands per circle."""
    if reference_id not in records:
        raise ValueError(f"Reference ID not found: {reference_id}")
    if kmer_size < 3 or any(len(value) < kmer_size for value in records.values()):
        raise ValueError("Anchor size must be at least three and fit every sequence")
    if any(not value or set(value) - IUPAC for value in records.values()):
        raise ValueError("Provide uppercase, ungapped IUPAC DNA sequences")
    ref = records[reference_id]
    reference_kmers = _kmers(ref, kmer_size)
    counts = Counter(reference_kmers + _kmers(reverse_complement(ref), kmer_size))
    candidates = {
        value for value in reference_kmers if counts[value] == 1 and not set(value) - set("ACGT")
    }
    for identifier, sequence in records.items():
        if identifier == reference_id:
            continue
        counts = Counter(_kmers(sequence, kmer_size) + _kmers(reverse_complement(sequence), kmer_size))
        candidates = {value for value in candidates if counts[value] == 1}
        if not candidates:
            raise ValueError(
                "No unique exact anchor shared by all circles. These inputs may be fragmented, "
                "repetitive or too divergent; use a different dataset or a shorter --anchor-size."
            )
    # Keep the deposited origin well inside the linear alignment rather than near its seam.
    positions = [index for index, value in enumerate(reference_kmers) if value in candidates]
    if not positions:
        raise ValueError("No unique exact A/C/G/T anchor: circular origin is ambiguous")
    return min(positions, key=lambda index: abs(index - len(ref) // 2))


def restore_reference_origin(
    aligned: Mapping[str, str], reference_id: str, reference_length: int, reference_rotation: int,
) -> dict[str, str]:
    """Rotate entire MSA columns; never infer an indel-sensitive query offset."""
    validate_alignment(aligned, reference_id)
    reference = aligned[reference_id]
    if len(reference.replace("-", "")) != reference_length:
        raise ValueError("Aligned reference length changed")
    target = (reference_length - reference_rotation) % reference_length
    base_index = 0
    cut = None
    for column, base in enumerate(reference):
        if base != "-":
            if base_index == target:
                cut = column
                break
            base_index += 1
    if cut is None:
        raise ValueError("Cannot locate reference origin in alignment")
    return {identifier: sequence[cut:] + sequence[:cut] for identifier, sequence in aligned.items()}


def align_circular(
    records: Mapping[str, str], reference_id: str, output: Path, mafft: str = "mafft",
    threads: int = 2, mafft_binaries: str | None = None, kmer_size: int = 31,
) -> dict[str, object]:
    """Align closely related complete circles; write auditable original coordinates.

    This is reference-guided circular normalization followed by MAFFT, not a
    globally optimal circular MSA or a reconstruction of ancestry.
    """
    if threads < 1:
        raise ValueError("threads must be positive")
    rotation = choose_common_anchor(records, reference_id, kmer_size)
    original_reference = records[reference_id]
    anchor_reference = original_reference[rotation:] + original_reference[:rotation]
    normalized: dict[str, str] = {}
    transformations: dict[str, object] = {}
    for identifier, sequence in records.items():
        result = normalize_to_reference(sequence, anchor_reference, kmer_size)
        normalized[identifier] = result.sequence
        transformation = asdict(result)
        del transformation["sequence"]
        transformations[identifier] = {
            **transformation, "input_length": len(sequence),
            "input_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        }
    output.mkdir(parents=True, exist_ok=True)
    write_fasta(normalized, output / "normalized.fasta")
    safe_ids = {f"S{index:06d}": identifier for index, identifier in enumerate(normalized)}
    command = [mafft, "--nuc", "--auto", "--thread", str(threads), "--inputorder"]
    environment = os.environ.copy()
    if mafft_binaries:
        environment["MAFFT_BINARIES"] = mafft_binaries
    with TemporaryDirectory(prefix="stump-of-life-") as temporary:
        temporary_path = Path(temporary)
        input_file = temporary_path / "input.fasta"
        output_file = temporary_path / "aligned.fasta"
        write_fasta({key: normalized[value] for key, value in safe_ids.items()}, input_file)
        with output_file.open("w") as stdout, (output / "mafft.log").open("w") as stderr:
            subprocess.run(
                [*command, str(input_file)], stdout=stdout, stderr=stderr,
                env=environment, check=True,
            )
        raw_alignment = read_fasta(output_file, aligned=True)
    if set(raw_alignment) != set(safe_ids):
        raise ValueError("MAFFT output does not contain exactly the input records")
    aligned = {safe_ids[key]: value for key, value in raw_alignment.items()}
    for identifier, sequence in aligned.items():
        if sequence.replace("-", "") != normalized[identifier]:
            raise ValueError(f"MAFFT changed input bases in {identifier}")
    aligned = restore_reference_origin(aligned, reference_id, len(original_reference), rotation)
    if aligned[reference_id].replace("-", "") != original_reference:
        raise ValueError("Reference origin restoration failed")
    write_fasta(aligned, output / "aligned.fasta")
    normalization = {
        "reference_id": reference_id, "reference_rotation_to_anchor": rotation,
        "anchor_sequence": anchor_reference[:kmer_size], "anchor_size": kmer_size,
        "rotation_convention": "Zero-based bases moved to the end after strand orientation; before MAFFT.",
        "ambiguity": "IUPAC retained; compatible ambiguity ignored in normalization edit distance, not called as SNP.",
        "samples": transformations,
    }
    (output / "normalization.json").write_text(json.dumps(normalization, indent=2) + "\n")
    summary = write_analysis(aligned, reference_id, output)
    summary.update({
        "method": "Reference-guided circular normalization followed by MAFFT",
        "mafft_command": command + ["<normalized input>"],
        "origin": "Original reference base 1; preceding insertions trail the last reference base.",
    })
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
