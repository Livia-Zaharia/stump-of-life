"""Render a gapped DNA alignment without confusing unknown calls with variants.

The SVG is a presentation drawing, not a cutting pattern or a CAD model. Full
mode preserves alignment-column angles, including reference-gap columns. Only
variable mode compresses columns. No sequence realignment happens here.
"""

from __future__ import annotations

import json
import math
from collections.abc import Collection, Mapping, Sequence
from pathlib import Path
from typing import Literal
from xml.etree import ElementTree as ET

_SVG = "http://www.w3.org/2000/svg"
_CALLED = frozenset("ACGT")
_AMBIGUOUS = frozenset("RYSWKMBDHVN")
_ALLOWED = _CALLED | _AMBIGUOUS | {"-"}
_COLORS = {"A": "#087f5b", "C": "#2563eb", "G": "#b76b00", "T": "#b52b59"}
_EVENTS = ("match", "snp", "deletion", "insertion", "unknown")
ET.register_namespace("", _SVG)


def _element(parent: ET.Element, tag: str, **attributes: str) -> ET.Element:
    return ET.SubElement(parent, f"{{{_SVG}}}{tag}", attributes)


def _text(
    parent: ET.Element,
    value: str,
    x: float,
    y: float,
    size: float = 16.0,
    anchor: str = "start",
    color: str = "#24323a",
) -> ET.Element:
    element = _element(
        parent,
        "text",
        x=f"{x:.3f}",
        y=f"{y:.3f}",
        fill=color,
        **{"font-size": f"{size:.3f}", "text-anchor": anchor},
    )
    element.text = value
    return element


def _is_xml_text(value: str) -> bool:
    return all(
        char in "\t\n\r"
        or 0x20 <= ord(char) <= 0xD7FF
        or 0xE000 <= ord(char) <= 0xFFFD
        or 0x10000 <= ord(char) <= 0x10FFFF
        for char in value
    )


def _event(reference: str, sample: str) -> str:
    """Call only differences supported by unambiguous bases or known gaps."""
    if reference in _AMBIGUOUS or sample in _AMBIGUOUS:
        return "unknown"
    if reference == sample:
        return "match"
    if reference == "-":
        return "insertion"
    if sample == "-":
        return "deletion"
    return "snp"


def _point(radius: float, angle: float) -> tuple[float, float]:
    return 500.0 + radius * math.sin(angle), 435.0 - radius * math.cos(angle)


def _draw_event(
    parent: ET.Element,
    event: str,
    base: str,
    radius: float,
    angle: float,
    half_size: float,
    attributes: dict[str, str],
    description: str,
) -> None:
    x, y = _point(radius, angle)
    if event in {"snp", "unknown"}:
        start_x, start_y = _point(radius - half_size, angle)
        end_x, end_y = _point(radius + half_size, angle)
        mark = _element(
            parent,
            "line",
            x1=f"{start_x:.4f}",
            y1=f"{start_y:.4f}",
            x2=f"{end_x:.4f}",
            y2=f"{end_y:.4f}",
            stroke="#8b929b" if event == "unknown" else _COLORS[base],
            **{"stroke-width": "1.35", "stroke-linecap": "round", **attributes},
        )
    elif event == "deletion":
        size = half_size * 0.75
        mark = _element(
            parent,
            "path",
            d=(
                f"M {x - size:.4f},{y - size:.4f} L {x + size:.4f},{y + size:.4f} "
                f"M {x - size:.4f},{y + size:.4f} L {x + size:.4f},{y - size:.4f}"
            ),
            stroke="#26343c",
            fill="none",
            **{"stroke-width": "1.15", "stroke-linecap": "round", **attributes},
        )
    else:
        size = half_size * 0.9
        mark = _element(
            parent,
            "polygon",
            points=(
                f"{x:.4f},{y - size:.4f} {x + size:.4f},{y:.4f} "
                f"{x:.4f},{y + size:.4f} {x - size:.4f},{y:.4f}"
            ),
            fill=_COLORS[base],
            stroke="#26343c",
            **{"stroke-width": "0.55", **attributes},
        )
    _element(mark, "title").text = description


def render_alignment_svg(
    records: Mapping[str, str],
    reference_id: str,
    output: Path,
    selected_ids: Sequence[str] | None = None,
    coordinate_mode: Literal["full", "variable"] = "full",
    diameter_mm: float = 500,
    excluded_columns: Collection[int] = (),
) -> dict[str, object]:
    """Write an SVG showing called changes against the innermost reference row.

    ``records`` must contain equally long, already gapped IUPAC DNA rows. The
    reference is always drawn, followed by selected rows in caller order; by
    default all other rows are selected. A matching canonical base is represented
    by its ring baseline. A/C/G/T changes are coloured ticks, deletions crosses,
    insertions diamonds, and comparisons involving ambiguous IUPAC symbols grey
    ticks. Ambiguous symbols never support a mutation call or establish a
    variable column by themselves.

    Column indices, including ``excluded_columns``, are zero based. Full mode
    uses each original alignment column divided by the original alignment length
    for its angle. Exclusions leave blank angular slots. Variable mode keeps
    columns containing at least two distinct known states (A/C/G/T or gap) among
    plotted rows, then spaces only those columns equally around the circle.

    ``diameter_mm`` sets the presentation canvas width, not a fabricated plate
    diameter. The canvas may extend vertically to accommodate the row legend.
    The returned JSON-serializable metadata is also embedded in the SVG. Event
    counts refer to rendered columns and are relative to the reference. Output
    is SVG only; optional raster conversion belongs to the caller.
    """
    if coordinate_mode not in {"full", "variable"}:
        raise ValueError("coordinate_mode must be 'full' or 'variable'.")
    if (
        not isinstance(diameter_mm, (int, float))
        or isinstance(diameter_mm, bool)
        or not math.isfinite(diameter_mm)
        or diameter_mm <= 0
    ):
        raise ValueError("diameter_mm must be a finite positive number.")
    if not records:
        raise ValueError("records must contain at least one aligned DNA sequence.")
    if reference_id not in records:
        raise ValueError(f"Unknown reference ID: {reference_id!r}.")
    sequences: dict[str, str] = {}
    for record_id, sequence in records.items():
        if not isinstance(record_id, str) or not record_id or not _is_xml_text(record_id):
            raise ValueError("Record IDs must be non-empty strings valid in XML text.")
        if not isinstance(sequence, str) or not sequence:
            raise ValueError(f"The aligned sequence for {record_id!r} is empty or not text.")
        normalized = sequence.upper()
        invalid = sorted(set(normalized) - _ALLOWED)
        if invalid:
            raise ValueError(
                f"The aligned sequence for {record_id!r} contains invalid DNA symbols: "
                f"{', '.join(repr(char) for char in invalid)}. Use IUPAC DNA and '-' gaps."
            )
        sequences[record_id] = normalized
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        detail = ", ".join(f"{key!r}: {len(value)}" for key, value in sequences.items())
        raise ValueError(f"Aligned rows must have equal lengths; received {detail}.")
    alignment_length = len(sequences[reference_id])
    if set(sequences[reference_id]) == {"-"}:
        raise ValueError("The reference row is entirely gaps; provide a reference containing DNA.")

    if selected_ids is None:
        chosen = list(records)
    else:
        if isinstance(selected_ids, (str, bytes)):
            raise ValueError("selected_ids must be a sequence of record IDs, not one string.")
        chosen = list(selected_ids)
        if any(not isinstance(record_id, str) for record_id in chosen):
            raise ValueError("Every selected ID must be a string.")
        if len(chosen) != len(set(chosen)):
            raise ValueError("selected_ids contains duplicate record IDs.")
        missing = [record_id for record_id in chosen if record_id not in records]
        if missing:
            raise ValueError(
                f"Unknown selected record IDs: {', '.join(repr(key) for key in missing)}."
            )
    plotted_ids = [reference_id] + [key for key in chosen if key != reference_id]
    excluded: set[int] = set()
    for column in excluded_columns:
        if not isinstance(column, int) or isinstance(column, bool):
            raise ValueError("excluded_columns must contain zero-based integer alignment indices.")
        if not 0 <= column < alignment_length:
            raise ValueError(
                f"Excluded column {column} is outside the alignment (0..{alignment_length - 1})."
            )
        excluded.add(column)
    available = [column for column in range(alignment_length) if column not in excluded]
    variable = [
        column
        for column in available
        if len({sequences[key][column] for key in plotted_ids} - _AMBIGUOUS) > 1
    ]
    columns = available if coordinate_mode == "full" else variable
    if not columns:
        if coordinate_mode == "variable":
            raise ValueError(
                "No informative variable columns remain among the selected rows. "
                "Unknown/ambiguous IUPAC calls alone are not mutations; try full mode "
                "or change the selection/exclusions."
            )
        raise ValueError("No alignment columns remain after applying excluded_columns.")

    reference = sequences[reference_id]
    reference_positions: list[int] = []
    reference_position = 0
    for base in reference:
        reference_position += base != "-"
        reference_positions.append(reference_position)
    counts = {key: {event: 0 for event in _EVENTS} for key in plotted_ids}
    for key in plotted_ids:
        for column in columns:
            counts[key][_event(reference[column], sequences[key][column])] += 1
    totals = {event: sum(counts[key][event] for key in plotted_ids[1:]) for event in _EVENTS}
    unknown_symbols = {
        key: sum(sequences[key][column] in _AMBIGUOUS for column in columns) for key in plotted_ids
    }
    metadata: dict[str, object] = {
        "reference_id": reference_id,
        "selected_ids": plotted_ids,
        "rows": len(plotted_ids),
        "columns": len(columns),
        "alignment_columns": alignment_length,
        "variable_columns": len(variable),
        "plotted_columns": columns,
        "excluded_columns": sorted(excluded),
        "mode": coordinate_mode,
        "coordinate_basis": "zero-based original multiple-alignment column",
        "angular_slots": alignment_length if coordinate_mode == "full" else len(columns),
        "diameter_mm": float(diameter_mm),
        "unknown_count": sum(unknown_symbols.values()),
        "unknown_symbols_by_row": unknown_symbols,
        "event_counts": counts,
        "comparison_totals": totals,
        "output": str(output),
    }
    mode_description = (
        "Full alignment mode: angles preserve original zero-based alignment-column positions, "
        "including insertions. Excluded columns retain blank angular slots."
        if coordinate_mode == "full"
        else "Compressed variable mode: only informative variable columns among the selected rows "
        "are spaced equally. Original genomic and alignment distances are not preserved."
    )
    description = (
        f"{mode_description} Reference {reference_id!r} is the innermost ring. "
        "Coloured A/C/G/T ticks are called substitutions; crosses are deletions; "
        "diamonds are insertions; grey ticks are unknown comparisons, not mutations. "
        "Matching called bases follow the ring baseline. This is a presentation SVG, not CAD."
    )
    legend_rows = math.ceil(len(plotted_ids) / 2)
    height = max(1000, 866 + legend_rows * 19 + 50)
    root = ET.Element(
        f"{{{_SVG}}}svg",
        {
            "width": f"{diameter_mm:g}mm",
            "height": f"{diameter_mm * height / 1000:g}mm",
            "viewBox": f"0 0 1000 {height}",
            "role": "img",
            "aria-labelledby": "alignment-title alignment-description",
            "font-family": "Arial, DejaVu Sans, sans-serif",
        },
    )
    _element(root, "title", id="alignment-title").text = "Stump of Life: circular DNA alignment"
    _element(root, "desc", id="alignment-description").text = description
    _element(root, "metadata").text = json.dumps(metadata, ensure_ascii=False, sort_keys=True)
    _element(root, "rect", width="1000", height=str(height), fill="#fcfbf7")
    _text(root, "STUMP OF LIFE", 500, 37, 24, "middle")
    mode_label = (
        "Full alignment coordinates" if coordinate_mode == "full" else "Compressed variable columns"
    )
    _text(
        root,
        f"{mode_label} | {len(plotted_ids)} rows | {len(columns):,} / {alignment_length:,} columns",
        500,
        63,
        15,
        "middle",
    )
    _text(
        root,
        f"{totals['snp']:,} substitutions   {totals['insertion']:,} insertions   "
        f"{totals['deletion']:,} deletions   {sum(unknown_symbols.values()):,} ambiguous calls",
        500,
        85,
        13,
        "middle",
        "#59656b",
    )
    rings = _element(root, "g", id="alignment-rings")
    spacing = 190.0 / max(1, len(plotted_ids) - 1)
    half_size = min(4.0, max(0.4, spacing * 0.30))
    for row_index, record_id in enumerate(plotted_ids):
        radius = 145.0 + row_index * spacing
        group = _element(
            rings,
            "g",
            id=f"ring-{row_index}",
            **{
                "data-record-id": record_id,
                "data-row": str(row_index),
                "data-radius": f"{radius:.4f}",
            },
        )
        _element(
            group,
            "circle",
            cx="500",
            cy="435",
            r=f"{radius:.4f}",
            fill="none",
            stroke="#667b80" if row_index == 0 else "#ccd3d2",
            **{"stroke-width": "1.3" if row_index == 0 else "0.65"},
        )
        _element(group, "title").text = f"Ring {row_index + 1}: {record_id}"
        for compressed_index, column in enumerate(columns):
            base = sequences[record_id][column]
            event = _event(reference[column], base)
            if event == "match":
                continue
            slot = column if coordinate_mode == "full" else compressed_index
            denominator = alignment_length if coordinate_mode == "full" else len(columns)
            angle = 2 * math.pi * slot / denominator
            attributes = {
                "data-event": event,
                "data-column": str(column),
                "data-base": base,
                "data-reference-base": reference[column],
                "data-reference-position": str(reference_positions[column]),
                "data-angle-deg": f"{360 * slot / denominator:.8f}",
            }
            coordinate = (
                f"insertion after reference position {reference_positions[column]}"
                if reference[column] == "-"
                else f"reference position {reference_positions[column]}"
            )
            _draw_event(
                group,
                event,
                base,
                radius,
                angle,
                half_size,
                attributes,
                f"{record_id}; alignment column {column} (zero based); {coordinate}; "
                f"{reference[column]} to {base}; {event}",
            )
    _text(root, "REFERENCE", 500, 416, 15, "middle", "#59656b")
    _text(root, reference_id, 500, 444, max(9.0, min(20.0, 240 / len(reference_id))), "middle")
    _text(root, "Ring 1; clockwise from the top", 500, 467, 12, "middle", "#59656b")
    _text(root, "A", 180, 804, 15, "middle", _COLORS["A"])
    _text(root, "C", 220, 804, 15, "middle", _COLORS["C"])
    _text(root, "G", 260, 804, 15, "middle", _COLORS["G"])
    _text(root, "T", 300, 804, 15, "middle", _COLORS["T"])
    _text(
        root,
        "tick = substitution   x = deletion   diamond = insertion   grey = unknown",
        330,
        804,
        12,
    )
    _text(root, "RINGS: INNER TO OUTER", 65, 844, 13, color="#59656b")
    for index, record_id in enumerate(plotted_ids):
        label_x = 65 if index % 2 == 0 else 525
        label_y = 868 + (index // 2) * 19
        label = _text(root, f"{index + 1}. {record_id}", label_x, label_y, 12)
        _element(label, "title").text = record_id
    _text(
        root,
        "Presentation drawing; not a fabrication pattern.",
        500,
        height - 20,
        12,
        "middle",
        "#59656b",
    )
    output = Path(output)
    if output.suffix.lower() != ".svg":
        raise ValueError(
            "render_alignment_svg writes SVG; use a .svg path and convert separately if needed."
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(root, space="  ")
    ET.ElementTree(root).write(output, encoding="utf-8", xml_declaration=True)
    return metadata
