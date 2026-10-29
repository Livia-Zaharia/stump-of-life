"""Export substitution-only circular perforation drawings in physical millimetres."""

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


def _element(tag: str, **attributes: str) -> ET.Element:
    return ET.Element(f"{{{_SVG}}}{tag}", attributes)


def _number(value: float) -> str:
    return f"{value:.9f}".rstrip("0").rstrip(".") or "0"


def _circle(parent: ET.Element, x: float, y: float, radius: float) -> ET.Element:
    return ET.SubElement(parent, f"{{{_SVG}}}circle", {
        "cx": _number(x), "cy": _number(y), "r": _number(radius),
    })


def _text(parent: ET.Element, x: float, y: float, value: str, size: float = 8) -> None:
    text = ET.SubElement(parent, f"{{{_SVG}}}text", {
        "x": _number(x), "y": _number(y), "font-size": _number(size),
        "font-family": "sans-serif", "fill": "#172a3a",
    })
    text.text = value


def _minimum_clearance(points: Sequence[tuple[float, float]], hole_diameter: float) -> float | None:
    if len(points) < 2:
        return None
    ordered = sorted(points)
    best_distance = math.inf
    # An x separation at least as large as the current best distance cannot
    # improve it. This is exact and prunes the common dense-circle case.
    for index, (x, y) in enumerate(ordered[:-1]):
        for other_index in range(index + 1, len(ordered)):
            other_x, other_y = ordered[other_index]
            if other_x - x >= best_distance:
                break
            best_distance = min(best_distance, math.hypot(other_x - x, other_y - y))
            if best_distance == 0:
                return -hole_diameter
    return best_distance - hole_diameter


def _write_svg(path: Path, root: ET.Element) -> None:
    ET.indent(root)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def export_perforation(
    records: Mapping[str, str],
    reference_id: str,
    output: Path,
    diameter_mm: float = 500,
    hole_diameter_mm: float = 5,
    min_web_mm: float = 1.5,
    coordinate_mode: Literal["full", "variable"] = "variable",
    selected_ids: Sequence[str] | None = None,
    excluded_columns: Collection[int] = (),
) -> dict[str, object]:
    """Write CAD and preview SVGs, hole provenance, and clearance measurements.

    ``excluded_columns`` uses zero-based MSA indices, matching
    ``informative_columns``. All exported column numbers are one-based.
    ``selected_ids`` may include the reference; it is always included once as
    the inner engraved guide. At least one non-reference row is required.
    Variable mode compresses invariant distances using the same informative
    columns across the reference and all selected rows, including indel columns.
    Only confidently called substitutions create holes in either mode.

    Clearance failures are exported with ``passes_requested_clearance=False``.
    These geometric measurements do not validate material strength or loading.
    """
    columns = validate_alignment(records, reference_id)
    for name, value in (
        ("diameter_mm", diameter_mm), ("hole_diameter_mm", hole_diameter_mm),
        ("min_web_mm", min_web_mm),
    ):
        if isinstance(value, bool) or not math.isfinite(value):
            raise ValueError(f"{name} must be finite.")
        if value < 0 or (name != "min_web_mm" and value == 0):
            raise ValueError(f"{name} must be positive (min_web_mm may be zero).")
    if hole_diameter_mm >= diameter_mm:
        raise ValueError("hole_diameter_mm must be smaller than diameter_mm.")
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
    modern = [identifier for identifier in selected if identifier != reference_id]
    if not modern:
        raise ValueError("Select at least one non-reference sample for the perforation drawing.")
    selected_records = {identifier: records[identifier] for identifier in [reference_id, *modern]}
    angular_columns = (
        list(range(columns)) if coordinate_mode == "full"
        else informative_columns(selected_records, exclusions)
    )
    description = (
        "Full MSA coordinates: every alignment column has an equal angular interval."
        if coordinate_mode == "full" else
        "Variable-column coordinates: selected informative MSA columns have equal angular "
        "intervals; invariant distances are compressed. Indel columns retain angular slots, "
        "but only confidently called substitutions create holes."
    )
    radius = diameter_mm / 2
    guide_radii = [0.55 * radius] + [
        radius * (0.83 if len(modern) == 1 else 0.70 + 0.26 * index / (len(modern) - 1))
        for index in range(len(modern))
    ]
    mapping = column_map(records[reference_id])
    holes: list[dict[str, object]] = []
    points: list[tuple[float, float]] = []
    for identifier, row_radius in zip(modern, guide_radii[1:], strict=True):
        for angular_index, column in enumerate(angular_columns):
            if column in exclusions:
                continue
            reference_base = records[reference_id][column]
            sample_base = records[identifier][column]
            event = call_state(reference_base, sample_base)
            if event != "substitution":
                continue
            angle_degrees = 360 * angular_index / len(angular_columns)
            angle = math.radians(angle_degrees - 90)
            x = radius + row_radius * math.cos(angle)
            y = radius + row_radius * math.sin(angle)
            points.append((x, y))
            holes.append({
                "sample_id": identifier, **mapping[column],
                "reference_base": reference_base, "sample_base": sample_base,
                "event": event, "angular_index": angular_index,
                "angle_degrees": angle_degrees, "radius_mm": row_radius,
                "x_mm": x, "y_mm": y, "hole_diameter_mm": hole_diameter_mm,
            })
    minimum_clearance = _minimum_clearance(points, hole_diameter_mm)
    edge_clearance = min(
        (radius - math.hypot(x - radius, y - radius) - hole_diameter_mm / 2 for x, y in points),
        default=None,
    )
    passes = all(
        clearance is None or clearance + 1e-9 >= min_web_mm
        for clearance in (minimum_clearance, edge_clearance)
    )
    output = output.resolve()
    comparison_counts = {
        identifier: dict(Counter(
            call_state(records[reference_id][column], records[identifier][column])
            for column in range(columns) if column not in exclusions
        ))
        for identifier in modern
    }
    geometry: dict[str, object] = {
        "design_status": "design_preview",
        "reference_id": reference_id, "selected_ids": modern,
        "coordinate_mode": coordinate_mode, "coordinate_description": description,
        "diameter_mm": diameter_mm, "hole_diameter_mm": hole_diameter_mm,
        "requested_min_web_mm": min_web_mm, "alignment_columns": columns,
        "angular_columns": len(angular_columns),
        "angular_alignment_columns": [column + 1 for column in angular_columns],
        "excluded_alignment_columns": sorted(column + 1 for column in exclusions),
        "hole_count": len(holes), "has_perforations": bool(holes),
        "comparison_counts_all_unexcluded_columns": comparison_counts,
        "min_hole_hole_clearance_mm": minimum_clearance,
        "edge_clearance_mm": edge_clearance,
        "passes_requested_clearance": passes,
        "validation_scope": (
            "Hole-edge to hole-edge and hole-edge to outer-edge distances only. "
            "No material, structural, sitting-load, or fabrication validation. "
            "Absent distances are null; a drawing with no holes passes spacing vacuously."
        ),
        "guides": [
            {"sample_id": identifier, "radius_mm": row_radius, "operation": "ENGRAVE"}
            for identifier, row_radius in zip([reference_id, *modern], guide_radii, strict=True)
        ],
        "files": {
            name: str(output / filename) for name, filename in (
                ("cut_svg", "cut.svg"), ("preview_svg", "preview.svg"),
                ("holes_csv", "holes.csv"), ("geometry_json", "geometry.json"),
            )
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    with (output / "holes.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "sample_id", *mapping[0], "reference_base", "sample_base", "event",
            "angular_index", "angle_degrees", "radius_mm", "x_mm", "y_mm", "hole_diameter_mm",
        ])
        writer.writeheader()
        writer.writerows(holes)
    (output / "geometry.json").write_text(
        json.dumps(geometry, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )

    drawing = _element(
        "svg", width=f"{_number(diameter_mm)}mm", height=f"{_number(diameter_mm)}mm",
        viewBox=f"0 0 {_number(diameter_mm)} {_number(diameter_mm)}",
    )
    ET.SubElement(drawing, f"{{{_SVG}}}title").text = "DNA perforation — design preview"
    ET.SubElement(drawing, f"{{{_SVG}}}desc").text = (
        f"Design preview, not fabrication approval. {description} "
        "Red CUT: outer boundary and substitution holes. Blue ENGRAVE: ring guides. "
        "No structural or sitting-load validation."
    )
    ET.SubElement(drawing, f"{{{_SVG}}}metadata").text = json.dumps({
        "design_status": "design_preview", "units": "mm",
        "passes_requested_clearance": passes, "coordinate_mode": coordinate_mode,
    })
    engraved = ET.SubElement(drawing, f"{{{_SVG}}}g", {
        "id": "ENGRAVE", "stroke": "#0000ff", "fill": "none", "stroke-width": "0.1",
    })
    for row_radius in guide_radii:
        _circle(engraved, radius, radius, row_radius)
    cuts = ET.SubElement(drawing, f"{{{_SVG}}}g", {
        "id": "CUT", "stroke": "#ff0000", "fill": "none", "stroke-width": "0.1",
    })
    _circle(cuts, radius, radius, radius)
    for x, y in points:
        _circle(cuts, x, y, hole_diameter_mm / 2)
    _write_svg(output / "cut.svg", drawing)

    legend_ids = [reference_id, *modern][:25]
    legend_width = max(170, max(len(identifier) for identifier in legend_ids) * 4.2 + 70)
    width = diameter_mm + legend_width + 80
    height = max(diameter_mm + 130, 130 + 11 * (len(legend_ids) + 1))
    preview = _element(
        "svg", width="1200", height=_number(1200 * height / width),
        viewBox=f"-30 -85 {_number(width)} {_number(height)}",
    )
    ET.SubElement(preview, f"{{{_SVG}}}rect", {
        "x": "-30", "y": "-85", "width": _number(width),
        "height": _number(height), "fill": "#fcfbf7",
    })
    _text(preview, 0, -55, "Circular DNA perforation — design preview", 12)
    _text(preview, 0, -38, f"{coordinate_mode.title()} coordinates | {len(holes)} substitution holes")
    _text(preview, 0, -26, (
        "Invariant distances compressed; equal angles for informative columns."
        if coordinate_mode == "variable" else "Equal angles for every MSA column."
    ))
    _text(preview, 0, -14, (
        "Requested geometric clearance: PASS" if passes else "Requested geometric clearance: FAIL"
    ))
    visible_guides = ET.SubElement(preview, f"{{{_SVG}}}g", {
        "stroke": "#6699dd", "fill": "none", "stroke-width": "0.45",
    })
    for row_radius in guide_radii:
        _circle(visible_guides, radius, radius, row_radius)
    _circle(visible_guides, radius, radius, radius).set("stroke", "#172a3a")
    visible_holes = ET.SubElement(preview, f"{{{_SVG}}}g", {"fill": "#172a3a"})
    for x, y in points:
        _circle(visible_holes, x, y, hole_diameter_mm / 2)
    legend_x = diameter_mm + 20
    _text(preview, legend_x, 8, "Engraved guides, inner to outer", 8.5)
    for index, identifier in enumerate(legend_ids):
        role = " (reference)" if index == 0 else ""
        _text(preview, legend_x, 22 + 11 * index, f"{index + 1}. {identifier}{role}")
    if len(modern) + 1 > len(legend_ids):
        _text(preview, legend_x, 22 + 11 * len(legend_ids), "Further guides: see geometry.json")
    _text(preview, 0, diameter_mm + 20, "Spacing check only; material and load capacity are unvalidated.")
    _write_svg(output / "preview.svg", preview)
    return geometry
