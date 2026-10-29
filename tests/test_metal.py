"""Perforations preserve substitution provenance and expose unsuitable spacing."""

import csv
import json
import math
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from stump_of_life.metal import export_perforation

SVG = {"svg": "http://www.w3.org/2000/svg"}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_full_coordinates_reveal_collisions_variable_coordinates_resolve_them(tmp_path: Path) -> None:
    records = {"reference": "A" * 1000, "sample": "CC" + "A" * 998}
    full = export_perforation(records, "reference", tmp_path / "full", coordinate_mode="full")
    variable = export_perforation(records, "reference", tmp_path / "variable")
    assert full["hole_count"] == variable["hole_count"] == 2
    assert full["passes_requested_clearance"] is False
    assert full["min_hole_hole_clearance_mm"] < 0
    assert variable["passes_requested_clearance"] is True
    assert variable["angular_alignment_columns"] == [1, 2]
    assert "compressed" in variable["coordinate_description"]
    rows = _rows(tmp_path / "full" / "holes.csv")
    assert [float(row["angle_degrees"]) for row in rows] == [0, 0.36]


def test_gaps_and_ambiguous_calls_are_never_holes(tmp_path: Path) -> None:
    records = {"reference": "A-CGN", "sample": "-TCNA"}
    result = export_perforation(records, "reference", tmp_path)
    assert result["hole_count"] == 0
    assert result["has_perforations"] is False
    assert result["min_hole_hole_clearance_mm"] is None
    assert result["edge_clearance_mm"] is None
    assert _rows(tmp_path / "holes.csv") == []
    json.dumps(result, allow_nan=False)


def test_one_modern_sample_has_a_finite_outer_ring_and_single_hole(tmp_path: Path) -> None:
    result = export_perforation({"reference": "ACGT", "one": "ATGT"}, "reference", tmp_path)
    assert result["hole_count"] == 1
    assert result["min_hole_hole_clearance_mm"] is None
    assert result["passes_requested_clearance"] is True
    assert result["guides"] == [
        {"sample_id": "reference", "radius_mm": 137.5, "operation": "ENGRAVE"},
        {"sample_id": "one", "radius_mm": 207.5, "operation": "ENGRAVE"},
    ]
    row = _rows(tmp_path / "holes.csv")[0]
    assert row["reference_position"] == "2"
    assert row["alignment_column"] == "2"
    assert row["reference_base"] == "C"
    assert row["sample_base"] == "T"
    assert row["event"] == "substitution"
    assert math.isfinite(float(row["x_mm"]))


def test_cut_svg_has_real_mm_dimensions_separate_operations_and_no_labels(tmp_path: Path) -> None:
    result = export_perforation({"ref": "ACGT", "sample<&": "ATGT"}, "ref", tmp_path)
    drawing = ET.parse(tmp_path / "cut.svg").getroot()
    assert drawing.attrib["width"] == "500mm"
    assert drawing.attrib["height"] == "500mm"
    assert drawing.attrib["viewBox"] == "0 0 500 500"
    assert drawing.find("svg:metadata", SVG) is not None
    assert "design preview" in drawing.find("svg:title", SVG).text.lower()
    assert drawing.findall(".//svg:text", SVG) == []
    groups = {group.attrib["id"]: group for group in drawing.findall("svg:g", SVG)}
    assert groups["CUT"].attrib["stroke"] == "#ff0000"
    assert groups["ENGRAVE"].attrib["stroke"] == "#0000ff"
    assert len(groups["CUT"]) == result["hole_count"] + 1
    assert len(groups["ENGRAVE"]) == 2
    assert json.loads((tmp_path / "geometry.json").read_text()) == result
    assert ET.parse(tmp_path / "preview.svg").findall(".//svg:text", SVG)


def test_selected_rows_share_one_angular_coordinate_system(tmp_path: Path) -> None:
    records = {"ref": "AAAA", "one": "CAAA", "two": "AAGA", "ignored": "AAAT"}
    result = export_perforation(records, "ref", tmp_path, selected_ids=["two", "one"])
    assert result["selected_ids"] == ["two", "one"]
    assert result["angular_alignment_columns"] == [1, 3]
    rows = _rows(tmp_path / "holes.csv")
    assert [(row["sample_id"], row["angle_degrees"]) for row in rows] == [
        ("two", "180.0"), ("one", "0.0"),
    ]
    assert result["guides"][0]["sample_id"] == "ref"


def test_exclusions_use_zero_based_indices_without_compressing_full_angles(tmp_path: Path) -> None:
    records = {"ref": "AAAA", "sample": "CCAA"}
    result = export_perforation(
        records, "ref", tmp_path, coordinate_mode="full", excluded_columns=[0],
    )
    assert result["hole_count"] == 1
    assert result["excluded_alignment_columns"] == [1]
    row = _rows(tmp_path / "holes.csv")[0]
    assert row["alignment_column"] == "2"
    assert float(row["angle_degrees"]) == 90


def test_negative_edge_clearance_is_exported_as_a_failure(tmp_path: Path) -> None:
    result = export_perforation(
        {"ref": "AAAA", "one": "CAAA", "two": "ACAA"}, "ref", tmp_path,
        hole_diameter_mm=30,
    )
    assert result["edge_clearance_mm"] < 0
    assert result["passes_requested_clearance"] is False
    assert (tmp_path / "cut.svg").exists()


@pytest.mark.parametrize("arguments", [
    {"diameter_mm": 0}, {"diameter_mm": -1}, {"diameter_mm": math.nan},
    {"hole_diameter_mm": 0}, {"hole_diameter_mm": math.inf},
    {"hole_diameter_mm": 500}, {"min_web_mm": -1}, {"min_web_mm": math.nan},
    {"coordinate_mode": "other"}, {"excluded_columns": [-1]},
    {"excluded_columns": [4]}, {"selected_ids": ["missing"]},
    {"selected_ids": ["sample", "sample"]}, {"selected_ids": []},
])
def test_invalid_geometry_and_selection_raise_before_export(tmp_path: Path, arguments: dict) -> None:
    with pytest.raises(ValueError):
        export_perforation({"ref": "AAAA", "sample": "CAAA"}, "ref", tmp_path, **arguments)
    assert not (tmp_path / "cut.svg").exists()


def test_invalid_alignment_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="equal lengths"):
        export_perforation({"ref": "AAAA", "sample": "AAA"}, "ref", tmp_path)
