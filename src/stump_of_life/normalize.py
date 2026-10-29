"""Conservatively orient circular, closely related DNA against a reference.

Rotation is the number of leading bases moved to the end *after* choosing the
reported strand. Input ambiguity symbols survive normalization and compatible
IUPAC bases (including an unknown ``N``) do not contribute to edit distance.
An insertion immediately before reference base one is placed at the end.
"""

from collections import Counter
from dataclasses import dataclass
from math import ceil
from typing import Any

import edlib

_BASES: dict[str, frozenset[str]] = {
    "A": frozenset("A"), "C": frozenset("C"),
    "G": frozenset("G"), "T": frozenset("T"),
    "R": frozenset("AG"), "Y": frozenset("CT"),
    "S": frozenset("CG"), "W": frozenset("AT"),
    "K": frozenset("GT"), "M": frozenset("AC"),
    "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_EQUALITIES: list[tuple[str, str]] = [
    (left, right)
    for left, bases in _BASES.items()
    for right, other_bases in _BASES.items()
    if left != right and bases & other_bases
]


@dataclass(frozen=True, slots=True)
class CircleNormalization:
    """The oriented sequence and its transformation from the cleaned input."""

    sequence: str
    strand: str
    rotation: int
    edit_distance: int
    anchor_count: int


def _clean_sequence(value: str, label: str) -> str:
    sequence = "".join(value.split()).upper()
    if not sequence:
        raise ValueError(f"The {label} sequence is empty.")
    invalid = sorted(set(sequence) - _BASES.keys())
    if invalid:
        raise ValueError(
            f"The {label} sequence contains non-IUPAC DNA symbols: "
            f"{', '.join(invalid)}. Provide unaligned DNA without gap characters."
        )
    return sequence


def _rotate(sequence: str, rotation: int) -> str:
    return sequence[rotation:] + sequence[:rotation]


def _unique_kmers(sequence: str, kmer_size: int) -> dict[str, int]:
    """Index exact, unambiguous k-mers occurring once on the circle."""
    extended = sequence + sequence[:kmer_size - 1]
    positions: dict[str, int] = {}
    repeated: set[str] = set()
    for position in range(len(sequence)):
        kmer = extended[position:position + kmer_size]
        if set(kmer) - {"A", "C", "G", "T"}:
            continue
        if kmer in positions:
            repeated.add(kmer)
        else:
            positions[kmer] = position
    return {kmer: position for kmer, position in positions.items() if kmer not in repeated}


def _align(query: str, target: str, mode: str, limit: int) -> dict[str, Any]:
    return edlib.align(
        query, target, mode=mode, task="locations", k=limit,
        additionalEqualities=_EQUALITIES,
    )


def normalize_to_reference(
    sequence: str, reference: str, kmer_size: int = 23,
) -> CircleNormalization:
    """Normalize a near-identical circular sequence without changing its bases.

    Unique A/C/G/T anchors determine plausible orientations. A shared unique
    reference-prefix anchor gives an exact origin immediately. Otherwise, a
    semiglobal alignment of the reference against two query copies finds origin
    candidates; full-circle global edit distances validate them. Up to five
    percent compatible edit distance (at least three edits for short inputs) is
    allowed. Ambiguous origins or strands raise ``ValueError``.

    This locates the origin; it does not construct a multiple alignment or
    restore deleted bases. A deletion crossing reference base one can make the
    origin intrinsically ambiguous. For variable reference beginnings, choosing
    a common exact anchor before calling this function gives reproducible
    multiple-alignment boundaries.
    """
    if not isinstance(kmer_size, int) or isinstance(kmer_size, bool) or kmer_size < 3:
        raise ValueError("kmer_size must be an integer of at least three.")
    cleaned = _clean_sequence(sequence, "query")
    ref = _clean_sequence(reference, "reference")
    if cleaned == ref:
        anchors = _unique_kmers(ref, kmer_size) if kmer_size <= len(ref) else {}
        return CircleNormalization(cleaned, "+", 0, 0, len(anchors))
    if kmer_size > min(len(cleaned), len(ref)):
        raise ValueError("The k-mer size exceeds a sequence length; use a smaller kmer_size.")

    limit = max(3, ceil(0.05 * max(len(cleaned), len(ref))))
    if abs(len(cleaned) - len(ref)) > limit:
        raise ValueError("The sequences are too dissimilar in length for circular normalization.")
    reference_kmers = _unique_kmers(ref, kmer_size)
    if not reference_kmers:
        raise ValueError(
            "The reference has no unique A/C/G/T anchors: its origin is ambiguous "
            "because it is repetitive or contains too many unknown bases."
        )

    results: list[CircleNormalization] = []
    prefix = ref[:kmer_size]
    prefix_is_unique = reference_kmers.get(prefix) == 0
    orientations = (("+", cleaned), ("-", cleaned.translate(_COMPLEMENT)[::-1]))
    for strand, oriented in orientations:
        query_kmers = _unique_kmers(oriented, kmer_size)
        shared = reference_kmers.keys() & query_kmers.keys()
        if len(shared) < 2:
            continue
        if prefix_is_unique and prefix in query_kmers:
            rotation = query_kmers[prefix]
            normalized = _rotate(oriented, rotation)
            distance = int(_align(ref, normalized, "NW", limit)["editDistance"])
            if distance >= 0:
                results.append(CircleNormalization(
                    normalized, strand, rotation, distance, len(shared),
                ))
            continue

        votes: Counter[int] = Counter(
            (query_kmers[kmer] - reference_kmers[kmer]) % len(oriented)
            for kmer in shared
        )
        circular_match = _align(ref, oriented + oriented, "HW", limit)
        if circular_match["editDistance"] < 0:
            continue
        locations = circular_match["locations"]
        semiglobal_origins = {
            int(start) % len(oriented)
            for start, end in locations
            if start is not None and end - start + 1 <= len(oriented)
        }
        # Origin sites in a repeated or unknown stretch can all match equally.
        if len(semiglobal_origins) > 32:
            raise ValueError(
                "The circular origin is ambiguous: many equally good alignments "
                "begin inside a repetitive or unknown region. Use a common unique anchor."
            )
        voted_origins = set(sorted(votes, key=lambda offset: (-votes[offset], offset))[:12])
        candidates = semiglobal_origins | voted_origins
        checked: list[CircleNormalization] = []
        for rotation in sorted(candidates):
            normalized = _rotate(oriented, rotation)
            distance = int(_align(ref, normalized, "NW", limit)["editDistance"])
            if distance >= 0:
                checked.append(CircleNormalization(
                    normalized, strand, rotation, distance, len(shared),
                ))
        # Upstream indels displace offset votes by at most the edit count. HW
        # alone may prefer an unphysical path consuming a query base twice.
        radius = min((item.edit_distance for item in checked), default=limit)
        nearby = {
            (offset + delta) % len(oriented)
            for offset in voted_origins
            for delta in range(-radius, radius + 1)
        } - candidates
        if len(nearby) > 2000:
            raise ValueError(
                "Too many plausible circular origins; choose a common unique "
                "reference anchor before normalizing this variable sequence."
            )
        for rotation in sorted(nearby):
            normalized = _rotate(oriented, rotation)
            distance = int(_align(ref, normalized, "NW", limit)["editDistance"])
            if distance >= 0:
                checked.append(CircleNormalization(
                    normalized, strand, rotation, distance, len(shared),
                ))
        if not checked:
            continue
        best_distance = min(item.edit_distance for item in checked)
        best = [item for item in checked if item.edit_distance == best_distance]
        if len(best) > 1:
            raise ValueError(
                "The circular origin is ambiguous: multiple rotations have the same "
                "validated edit distance. Use a common unique reference anchor."
            )
        results.extend(best)

    if not results:
        raise ValueError(
            "No reliable circular normalization: sequences are too dissimilar "
            "(more than five percent edits), repetitive, or share too few unique "
            "A/C/G/T anchors. Try a smaller kmer_size only for short sequences."
        )
    best_distance = min(item.edit_distance for item in results)
    best_results = [item for item in results if item.edit_distance == best_distance]
    if len(best_results) != 1:
        raise ValueError(
            "The strand is ambiguous: forward and reverse-complement orientations "
            "have the same validated edit distance."
        )
    return best_results[0]
