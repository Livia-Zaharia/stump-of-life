"""Small, strict FASTA utilities with no implicit biological transformations."""

from collections.abc import Mapping
from pathlib import Path

IUPAC = frozenset("ACGTRYSWKMBDHVN")


def read_fasta(path: Path, *, aligned: bool = False) -> dict[str, str]:
    """Read unique first-word FASTA IDs, preserving order and IUPAC symbols."""
    records: dict[str, str] = {}
    identifier: str | None = None
    parts: list[str] = []

    def save() -> None:
        if identifier is None:
            return
        sequence = "".join(parts).upper()
        if not sequence:
            raise ValueError(f"Empty FASTA record: {identifier}")
        invalid = set(sequence) - IUPAC - ({"-"} if aligned else set())
        if invalid:
            raise ValueError(f"Invalid DNA symbols in {identifier}: {sorted(invalid)}")
        records[identifier] = sequence

    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            save()
            words = line[1:].split()
            if not words:
                raise ValueError(f"Missing FASTA ID on line {line_number}")
            identifier = words[0]
            if identifier in records:
                raise ValueError(f"Duplicate FASTA ID: {identifier}")
            parts = []
        elif identifier is None:
            raise ValueError(f"Sequence before a FASTA header on line {line_number}")
        else:
            parts.append("".join(line.split()))
    save()
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    if aligned and len({len(sequence) for sequence in records.values()}) != 1:
        raise ValueError("Aligned FASTA rows must have equal lengths")
    return records


def write_fasta(records: Mapping[str, str], path: Path) -> None:
    """Write records without silently truncating identifiers."""
    lines: list[str] = []
    for identifier, sequence in records.items():
        if not identifier or any(character.isspace() for character in identifier):
            raise ValueError("FASTA IDs must be nonempty and contain no whitespace")
        lines.append(f">{identifier}\n")
        lines.extend(f"{sequence[index:index + 80]}\n" for index in range(0, len(sequence), 80))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines))
