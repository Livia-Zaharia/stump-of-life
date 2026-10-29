"""Export aligned substitutions as surface marks on a solid circular cap."""

import csv
import json
import math
from collections import Counter
from collections.abc import Collection, Mapping, Sequence
from pathlib import Path
from typing import Literal
from xml.etree import ElementTree as ET

from stump_of_life.analysis import (
    call_state,
    column_map,
    informative_columns,
    validate_alignment,
)

_SVG = "http://www.w3.org/2000/svg"
ET.register_namespace("", _SVG)


def _number(value: float) -> str:
    return f"{value:.9f}".rstrip("0").rstrip(".") or "0"


def _circle(
    parent: ET.Element, x: float, y: float, radius: float,
    attributes: Mapping[str, str] | None = None,
) -> ET.Element:
    return ET.SubElement(parent, f"{{{_SVG}}}circle", {
        "cx": _number(x), "cy": _number(y), "r": _number(radius), **(attributes or {}),
    })


def _minimum_mark_gap(
    points: Sequence[tuple[float, float]], mark_diameter: float,
) -> float | None:
    if len(points) < 2:
        return None
    ordered = sorted(points)
    nearest = math.inf
    for index, (x, y) in enumerate(ordered[:-1]):
        for other_index in range(index + 1, len(ordered)):
            other_x, other_y = ordered[other_index]
            if other_x - x >= nearest:
                break
            nearest = min(nearest, math.hypot(other_x - x, other_y - y))
            if nearest == 0:
                return -mark_diameter
    return nearest - mark_diameter


def export_engraving(
    records: Mapping[str, str],
    reference_id: str,
    output: Path,
    diameter_mm: float = 500,
    mark_diameter_mm: float = 1.2,
    min_gap_mm: float = 0.2,
    coordinate_mode: Literal["full", "variable"] = "variable",
    selected_ids: Sequence[str] | None = None,
    excluded_columns: Collection[int] = (),
) -> dict[str, object]:
    """Export a solid cap, substitution marks, their provenance and spacing.

    The single CUT loop is the outer boundary. All row guides and filled dots
    belong to ENGRAVE; dots never become cutting loops or holes. The reference
    guide is at 0.55R; other rows occupy 0.70R–0.96R (0.83R for a single row).
    These are generic drawing radii, without an implied date or ancestry order.

    Input exclusions are zero-based MSA indices, matching ``informative_columns``;
    exported columns are one-based. Full mode keeps every column's angular slot,
    including excluded columns. Variable mode shares informative columns across
    the reference and selected rows, compressing invariant distances. Gaps and
    unknowns are counted but produce no substitution marks.

    Spacing failures are exported for review with the pass flag set to False.
    Distances measure dot-edge separation, not remaining structural material.
    """
    columns = validate_alignment(records, reference_id)
    for name, value in (
        ("diameter_mm", diameter_mm), ("mark_diameter_mm", mark_diameter_mm),
        ("min_gap_mm", min_gap_mm),
    ):
        if isinstance(value, bool) or not math.isfinite(value):
            raise ValueError(f"{name} must be finite.")
        if value < 0 or (name != "min_gap_mm" and value == 0):
            raise ValueError(f"{name} must be positive (min_gap_mm may be zero).")
    if mark_diameter_mm >= diameter_mm:
        raise ValueError("mark_diameter_mm must be smaller than diameter_mm.")
    if coordinate_mode not in {"full", "variable"}:
        raise ValueError("coordinate_mode must be 'full' or 'variable'.")
    exclusions = set(excluded_columns)
    if any(
        not isinstance(column, int) or isinstance(column, bool) or not 0 <= column < columns
        for column in exclusions
    ):
        raise ValueError("excluded_columns must contain zero-based alignment indices in range.")
    selected = list(records) if selected_ids is None else list(selected_ids)
    if len(selected) != len(set(selected)):
        raise ValueError("selected_ids contains duplicate sample identifiers.")
    missing = [identifier for identifier in selected if identifier not in records]
    if missing:
        raise ValueError(f"Selected sample identifiers not found: {', '.join(missing)}")
    samples = [identifier for identifier in selected if identifier != reference_id]
    if not samples:
        raise ValueError("Select at least one non-reference sample for the engraving drawing.")
    selected_records = {identifier: records[identifier] for identifier in [reference_id, *samples]}
    angular_columns = (
        list(range(columns)) if coordinate_mode == "full"
        else informative_columns(selected_records, exclusions)
    )
    coordinate_description = (
        "Full MSA coordinates: every column has an equal angular interval; excluded "
        "columns keep blank slots."
        if coordinate_mode == "full" else
        "Variable-column coordinates: selected informative columns share equal angular "
        "intervals across rows; invariant distances are compressed. Indel columns retain "
        "slots, but only confidently called substitutions receive engraved marks."
    )
    radius = diameter_mm / 2
    radii = [0.55 * radius] + [
        radius * (0.83 if len(samples) == 1 else 0.70 + 0.26 * index / (len(samples) - 1))
        for index in range(len(samples))
    ]
    mapping = column_map(records[reference_id])
    marks: list[dict[str, object]] = []
    points: list[tuple[float, float]] = []
    for identifier, row_radius in zip(samples, radii[1:], strict=True):
        for angular_index, column in enumerate(angular_columns):
            if column in exclusions:
                continue
            reference_base = records[reference_id][column]
            sample_base = records[identifier][column]
            if call_state(reference_base, sample_base) != "substitution":
                continue
            degrees = 360 * angular_index / len(angular_columns)
            angle = math.radians(degrees - 90)
            x = radius + row_radius * math.cos(angle)
            y = radius + row_radius * math.sin(angle)
            points.append((x, y))
            marks.append({
                "sample_id": identifier, **mapping[column],
                "reference_base": reference_base, "sample_base": sample_base,
                "event": "substitution", "operation": "ENGRAVE",
                "angular_index": angular_index, "angle_degrees": degrees,
                "radius_mm": row_radius, "x_mm": x, "y_mm": y,
                "mark_diameter_mm": mark_diameter_mm,
            })
    comparison_counts = {
        identifier: dict(Counter(
            call_state(records[reference_id][column], records[identifier][column])
            for column in range(columns) if column not in exclusions
        ))
        for identifier in samples
    }
    totals: Counter[str] = Counter()
    for counts in comparison_counts.values():
        totals.update(counts)
    mark_gap = _minimum_mark_gap(points, mark_diameter_mm)
    edge_clearance = min(
        (radius - math.hypot(x - radius, y - radius) - mark_diameter_mm / 2 for x, y in points),
        default=None,
    )
    passes = all(
        clearance is None or clearance + 1e-9 >= min_gap_mm
        for clearance in (mark_gap, edge_clearance)
    )
    output = output.resolve()
    geometry: dict[str, object] = {
        "design_status": "design_preview", "geometry_kind": "surface_engraving",
        "reference_id": reference_id, "selected_ids": samples,
        "diameter_mm": diameter_mm, "mark_diameter_mm": mark_diameter_mm,
        "requested_min_gap_mm": min_gap_mm, "coordinate_mode": coordinate_mode,
        "coordinate_description": coordinate_description, "alignment_columns": columns,
        "angular_columns": len(angular_columns),
        "angular_alignment_columns": [column + 1 for column in angular_columns],
        "excluded_alignment_columns": sorted(column + 1 for column in exclusions),
        "mark_count": len(marks), "has_marks": bool(marks),
        "surface_has_perforations": False, "cut_boundary_count": 1,
        "comparison_counts_all_unexcluded_columns": comparison_counts,
        "omitted_event_counts": {state: totals[state] for state in ("insertion", "deletion", "unknown")},
        "min_mark_mark_gap_mm": mark_gap, "edge_clearance_mm": edge_clearance,
        "passes_requested_clearance": passes,
        "validation_scope": (
            "Geometric separation between filled-dot edges and from dot edges to the "
            "outer boundary only; intentional intersections with engraved guides are excluded. "
            "These are readability/density measurements, not structural webs, material, load "
            "or manufacturing validation. Null distances mean no applicable pair or mark. "
            "A drawing without marks passes spacing vacuously."
        ),
        "guides": [
            {"sample_id": identifier, "radius_mm": row_radius, "operation": "ENGRAVE"}
            for identifier, row_radius in zip([reference_id, *samples], radii, strict=True)
        ],
        "files": {name: str(output / filename) for name, filename in (
            ("engraving_svg", "engraving.svg"), ("marks_csv", "marks.csv"),
            ("geometry_json", "geometry.json"),
        )},
    }
    output.mkdir(parents=True, exist_ok=True)
    with (output / "marks.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "sample_id", *mapping[0], "reference_base", "sample_base", "event", "operation",
            "angular_index", "angle_degrees", "radius_mm", "x_mm", "y_mm", "mark_diameter_mm",
        ])
        writer.writeheader()
        writer.writerows(marks)
    (output / "geometry.json").write_text(
        json.dumps(geometry, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    drawing = ET.Element(f"{{{_SVG}}}svg", {
        "width": f"{_number(diameter_mm)}mm", "height": f"{_number(diameter_mm)}mm",
        "viewBox": f"0 0 {_number(diameter_mm)} {_number(diameter_mm)}",
    })
    ET.SubElement(drawing, f"{{{_SVG}}}title").text = "Solid DNA cap engraving — design preview"
    ET.SubElement(drawing, f"{{{_SVG}}}desc").text = (
        f"{coordinate_description} CUT contains the outer perimeter only. ENGRAVE contains "
        "closed guide rings and filled substitution dots; no dot is a through-hole or cut loop. "
        "Dimensions are millimetres. Mark spacing is a density/readability check, not "
        "structural or fabrication approval."
    )
    ET.SubElement(drawing, f"{{{_SVG}}}metadata").text = json.dumps(geometry, allow_nan=False)
    cuts = ET.SubElement(drawing, f"{{{_SVG}}}g", {
        "id": "CUT", "stroke": "#ff0000", "stroke-width": "0.1", "fill": "none",
    })
    _circle(cuts, radius, radius, radius, {"data-operation": "CUT"})
    engraved = ET.SubElement(drawing, f"{{{_SVG}}}g", {"id": "ENGRAVE"})
    guides = ET.SubElement(engraved, f"{{{_SVG}}}g", {
        "id": "GUIDES", "stroke": "#0000ff", "stroke-width": "0.1", "fill": "none",
    })
    for identifier, row_radius in zip([reference_id, *samples], radii, strict=True):
        _circle(guides, radius, radius, row_radius, {
            "data-operation": "ENGRAVE", "data-sample-id": identifier,
        })
    dots = ET.SubElement(engraved, f"{{{_SVG}}}g", {
        "id": "SUBSTITUTIONS", "fill": "#172a3a", "stroke": "none",
    })
    for mark in marks:
        _circle(dots, float(mark["x_mm"]), float(mark["y_mm"]), mark_diameter_mm / 2, {
            "data-operation": "ENGRAVE", "data-sample-id": str(mark["sample_id"]),
            "data-alignment-column": str(mark["alignment_column"]),
            "data-reference-position": str(mark["reference_position"]),
            "data-reference-base": str(mark["reference_base"]),
            "data-sample-base": str(mark["sample_base"]),
        })
    ET.indent(drawing)
    ET.ElementTree(drawing).write(
        output / "engraving.svg", encoding="utf-8", xml_declaration=True,
    )
    return geometry
