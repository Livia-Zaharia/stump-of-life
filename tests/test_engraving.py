"""Surface dots must remain distinct from through-cut fabrication paths."""

import csv
import json
import math
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from stump_of_life.engraving import export_engraving

SVG = {"svg": "http://www.w3.org/2000/svg"}


def _marks(output: Path) -> list[dict[str, str]]:
    with (output / "marks.csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_cut_contains_only_outer_boundary_and_dots_are_surface_engraving(tmp_path: Path) -> None:
    result = export_engraving({"ref": "AAAA", "sample": "CAGA"}, "ref", tmp_path)
    drawing = ET.parse(tmp_path / "engraving.svg").getroot()
    assert drawing.attrib["width"] == drawing.attrib["height"] == "500mm"
    assert drawing.attrib["viewBox"] == "0 0 500 500"
    groups = {group.attrib["id"]: group for group in drawing.findall("svg:g", SVG)}
    cuts = groups["CUT"].findall(".//svg:circle", SVG)
    assert len(cuts) == 1
    assert cuts[0].attrib["r"] == "250"
    assert len(list(groups["CUT"])) == 1
    assert groups["CUT"].attrib["stroke"] == "#ff0000"
    dots = groups["ENGRAVE"].findall("svg:g[@id='SUBSTITUTIONS']/svg:circle", SVG)
    assert len(dots) == result["mark_count"] == 2
    assert all(dot.attrib["r"] == "0.6" for dot in dots)
    assert all(dot.attrib["data-operation"] == "ENGRAVE" for dot in dots)
    assert groups["ENGRAVE"].find("svg:g[@id='SUBSTITUTIONS']", SVG).attrib["stroke"] == "none"
    assert groups["ENGRAVE"].find("svg:g[@id='SUBSTITUTIONS']", SVG).attrib["fill"] != "none"
    assert len(groups["ENGRAVE"].findall("svg:g[@id='GUIDES']/svg:circle", SVG)) == 2
    assert result["surface_has_perforations"] is False
    assert result["geometry_kind"] == "surface_engraving"
    assert {path.name for path in tmp_path.iterdir()} == {
        "engraving.svg", "marks.csv", "geometry.json",
    }
    assert json.loads((tmp_path / "geometry.json").read_text()) == result


def test_unknowns_and_indels_are_counted_without_becoming_dots(tmp_path: Path) -> None:
    result = export_engraving({"ref": "A-CGNAT", "sample": "-TCTANC"}, "ref", tmp_path)
    assert result["mark_count"] == 2
    assert result["omitted_event_counts"] == {"insertion": 1, "deletion": 1, "unknown": 2}
    assert result["comparison_counts_all_unexcluded_columns"]["sample"] == {
        "deletion": 1, "insertion": 1, "match": 1, "substitution": 2, "unknown": 2,
    }
    rows = _marks(tmp_path)
    assert [(row["alignment_column"], row["reference_position"]) for row in rows] == [
        ("4", "3"), ("7", "6"),
    ]
    assert [(row["reference_base"], row["sample_base"]) for row in rows] == [("G", "T"), ("T", "C")]
    assert all(row["event"] == "substitution" and row["operation"] == "ENGRAVE" for row in rows)


def test_all_unknowns_give_no_dots_and_json_safe_null_distances(tmp_path: Path) -> None:
    result = export_engraving({"ref": "ACGT", "sample": "NNNN"}, "ref", tmp_path)
    assert result["mark_count"] == result["angular_columns"] == 0
    assert result["has_marks"] is False
    assert result["omitted_event_counts"]["unknown"] == 4
    assert result["min_mark_mark_gap_mm"] is None
    assert result["edge_clearance_mm"] is None
    assert _marks(tmp_path) == []
    json.dumps(result, allow_nan=False)


def test_full_coordinate_clusters_fail_while_variable_coordinates_separate(tmp_path: Path) -> None:
    records = {"ref": "A" * 10000, "sample": "CC" + "A" * 9998}
    full = export_engraving(records, "ref", tmp_path / "full", coordinate_mode="full")
    variable = export_engraving(records, "ref", tmp_path / "variable")
    assert full["mark_count"] == variable["mark_count"] == 2
    assert full["passes_requested_clearance"] is False
    assert full["min_mark_mark_gap_mm"] < 0
    assert variable["passes_requested_clearance"] is True
    assert variable["angular_alignment_columns"] == [1, 2]
    assert "compressed" in variable["coordinate_description"]
    assert (tmp_path / "full" / "engraving.svg").exists()


def test_dot_spacing_includes_neighbouring_tracks(tmp_path: Path) -> None:
    records = {"ref": "AA", **{f"sample_{index}": "CA" for index in range(100)}}
    result = export_engraving(records, "ref", tmp_path)
    assert result["mark_count"] == 100
    assert result["min_mark_mark_gap_mm"] == pytest.approx(65 / 99 - 1.2)
    assert result["passes_requested_clearance"] is False


def test_single_sample_and_selection_preserve_reference_and_column_identity(tmp_path: Path) -> None:
    result = export_engraving(
        {"ref": "AAAA", "selected": "AACA", "ignored": "CAAA"}, "ref", tmp_path,
        selected_ids=["selected"],
    )
    assert result["selected_ids"] == ["selected"]
    assert result["angular_alignment_columns"] == [3]
    assert [guide["radius_mm"] for guide in result["guides"]] == [137.5, 207.5]
    assert result["min_mark_mark_gap_mm"] is None
    assert result["passes_requested_clearance"] is True
    row = _marks(tmp_path)[0]
    assert (row["sample_id"], row["alignment_column"], row["reference_position"]) == (
        "selected", "3", "3",
    )
    assert float(row["angle_degrees"]) == 0
    assert math.isfinite(float(row["y_mm"]))


def test_excluded_full_columns_leave_their_original_angular_slots(tmp_path: Path) -> None:
    result = export_engraving(
        {"ref": "AAAA", "sample": "CCAA"}, "ref", tmp_path,
        coordinate_mode="full", excluded_columns=[0],
    )
    assert result["angular_alignment_columns"] == [1, 2, 3, 4]
    assert result["excluded_alignment_columns"] == [1]
    assert result["mark_count"] == 1
    row = _marks(tmp_path)[0]
    assert row["alignment_column"] == "2"
    assert float(row["angle_degrees"]) == 90
    assert result["comparison_counts_all_unexcluded_columns"]["sample"]["substitution"] == 1


def test_sample_identifiers_are_xml_escaped_without_becoming_machine_ids(tmp_path: Path) -> None:
    identifier = 'sample <& "foreign">'
    export_engraving({"ref": "AA", identifier: "CA"}, "ref", tmp_path)
    drawing = ET.parse(tmp_path / "engraving.svg").getroot()
    dot = drawing.find("svg:g[@id='ENGRAVE']/svg:g[@id='SUBSTITUTIONS']/svg:circle", SVG)
    assert dot.attrib["data-sample-id"] == identifier
    assert _marks(tmp_path)[0]["sample_id"] == identifier
    assert {node.attrib["id"] for node in drawing.findall(".//svg:g", SVG)} == {
        "CUT", "ENGRAVE", "GUIDES", "SUBSTITUTIONS",
    }


def test_marks_outside_perimeter_export_with_a_failed_gap_report(tmp_path: Path) -> None:
    result = export_engraving(
        {"ref": "AA", "one": "CA", "two": "AC"}, "ref", tmp_path,
        diameter_mm=20, mark_diameter_mm=2,
    )
    assert result["edge_clearance_mm"] < 0
    assert result["passes_requested_clearance"] is False
    assert (tmp_path / "engraving.svg").exists()


@pytest.mark.parametrize("arguments", [
    {"diameter_mm": 0}, {"diameter_mm": -1}, {"diameter_mm": math.nan},
    {"mark_diameter_mm": 0}, {"mark_diameter_mm": math.inf},
    {"mark_diameter_mm": 500}, {"min_gap_mm": -1}, {"min_gap_mm": math.nan},
    {"coordinate_mode": "other"}, {"excluded_columns": [-1]},
    {"excluded_columns": [4]}, {"excluded_columns": [True]},
    {"selected_ids": ["missing"]}, {"selected_ids": []},
    {"selected_ids": ["sample", "sample"]},
])
def test_invalid_arguments_raise_before_writing(tmp_path: Path, arguments: dict) -> None:
    with pytest.raises(ValueError):
        export_engraving({"ref": "AAAA", "sample": "CAAA"}, "ref", tmp_path, **arguments)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("records", [
    {"ref": "", "sample": ""},
    {"ref": "----", "sample": "ACGT"},
    {"ref": "AAAA", "sample": "AAA"},
])
def test_empty_or_invalid_reference_alignments_are_rejected(tmp_path: Path, records: dict) -> None:
    with pytest.raises(ValueError):
        export_engraving(records, "ref", tmp_path)
