"""Keep alignment columns, reference coordinates and uncertain calls distinct."""

import csv
from collections import Counter
from collections.abc import Collection, Mapping
from pathlib import Path

KNOWN = frozenset("ACGT")


def call_state(reference: str, query: str) -> str:
    """Classify a column; ambiguous calls never become substitutions."""
    if reference not in KNOWN | {"-"} or query not in KNOWN | {"-"}:
        return "unknown"
    if reference == query:
        return "match"
    if reference == "-":
        return "insertion"
    if query == "-":
        return "deletion"
    return "substitution"


def validate_alignment(records: Mapping[str, str], reference_id: str) -> int:
    if reference_id not in records:
        raise ValueError(f"Reference ID not found: {reference_id}")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1 or not next(iter(lengths)):
        raise ValueError("Alignment rows must be nonempty and have equal lengths")
    if not records[reference_id].replace("-", ""):
        raise ValueError("Reference must contain DNA bases, not only gaps")
    invalid = set("".join(records.values())) - set("ACGTRYSWKMBDHVN-")
    if invalid:
        raise ValueError(f"Invalid aligned DNA symbols: {sorted(invalid)}")
    return lengths.pop()


def column_map(reference: str) -> list[dict[str, int | None]]:
    """One-based columns and positions; insertion columns carry flank and rank.

    In imported alignments with leading gaps, their seam ranks follow the
    trailing gaps so every circular insertion address is unique.
    """
    result: list[dict[str, int | None]] = []
    position = 0
    insertion = 0
    reference_length = len(reference.replace("-", ""))
    if not reference_length:
        raise ValueError("Reference must contain DNA bases, not only gaps")
    trailing_gaps = len(reference) - len(reference.rstrip("-"))
    for column, base in enumerate(reference, 1):
        if base == "-":
            insertion += 1
            result.append({
                "alignment_column": column,
                "reference_position": None,
                "insertion_after_reference_position": position or reference_length,
                "insertion_rank": insertion + (trailing_gaps if position == 0 else 0),
            })
        else:
            position += 1
            insertion = 0
            result.append({
                "alignment_column": column,
                "reference_position": position,
                "insertion_after_reference_position": None,
                "insertion_rank": None,
            })
    return result


def informative_columns(
    records: Mapping[str, str], excluded: Collection[int] = (),
) -> list[int]:
    """Columns with two distinct confidently known states (including gaps)."""
    exclusions = set(excluded)
    return [
        column for column, bases in enumerate(zip(*records.values(), strict=True))
        if column not in exclusions and len(set(bases) & (KNOWN | {"-"})) > 1
    ]


def write_analysis(records: Mapping[str, str], reference_id: str, output: Path) -> dict[str, object]:
    columns = validate_alignment(records, reference_id)
    mapping = column_map(records[reference_id])
    output.mkdir(parents=True, exist_ok=True)
    with (output / "columns.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(mapping[0]))
        writer.writeheader()
        writer.writerows(mapping)
    counts: dict[str, dict[str, int]] = {}
    with (output / "variants.csv").open("w", newline="") as handle:
        fields = ["sample_id", *mapping[0], "reference_base", "sample_base", "state"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for identifier, sequence in records.items():
            row_counts: Counter[str] = Counter()
            for index, (ref, query) in enumerate(zip(records[reference_id], sequence, strict=True)):
                state = call_state(ref, query)
                row_counts[state] += 1
                if state != "match":
                    writer.writerow({
                        "sample_id": identifier, **mapping[index],
                        "reference_base": ref, "sample_base": query, "state": state,
                    })
            counts[identifier] = dict(row_counts)
    return {
        "reference_id": reference_id,
        "reference_length": len(records[reference_id].replace("-", "")),
        "alignment_columns": columns,
        "samples": len(records),
        "variable_columns": len(informative_columns(records)),
        "counts_by_sample": counts,
        "coordinates": "One-based columns and positions in the selected reference; not ancestral states.",
    }
