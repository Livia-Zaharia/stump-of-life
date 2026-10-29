"""Scientific coordinate and calling invariants of the SVG renderer."""

import json
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from stump_of_life.render import render_alignment_svg

SVG = "{http://www.w3.org/2000/svg}"


def events(svg_path: Path, record_id: str) -> list[ET.Element]:
    root = ET.parse(svg_path).getroot()
    group = next(item for item in root.iter(f"{SVG}g") if item.get("data-record-id") == record_id)
    return [item for item in group if item.get("data-event") is not None]


def test_unknown_iupac_calls_are_grey_and_never_substitutions(tmp_path: Path) -> None:
    output = tmp_path / "unknown.svg"
    ambiguous = "RYSWKMBDHVN"
    metadata = render_alignment_svg(
        {"ref": "A" * len(ambiguous), "sample": ambiguous},
        "ref",
        output,
    )
    marks = events(output, "sample")
    assert len(marks) == len(ambiguous)
    assert {mark.get("data-event") for mark in marks} == {"unknown"}
    assert {mark.get("stroke") for mark in marks} == {"#8b929b"}
    assert metadata["unknown_count"] == len(ambiguous)
    assert metadata["comparison_totals"]["snp"] == 0
    with pytest.raises(ValueError, match="No informative variable columns"):
        render_alignment_svg(
            {"ref": "A" * len(ambiguous), "sample": ambiguous},
            "ref",
            tmp_path / "unknown-variable.svg",
            coordinate_mode="variable",
        )


def test_unknown_reference_does_not_turn_called_sample_into_snp(tmp_path: Path) -> None:
    output = tmp_path / "unknown-reference.svg"
    metadata = render_alignment_svg({"ref": "AN", "sample": "AC"}, "ref", output)
    assert metadata["event_counts"]["sample"]["snp"] == 0
    assert metadata["event_counts"]["sample"]["unknown"] == 1
    assert metadata["unknown_count"] == 1
    assert events(output, "sample")[0].get("data-event") == "unknown"


def test_substitutions_are_ticks_coloured_by_called_sample_base(tmp_path: Path) -> None:
    output = tmp_path / "called.svg"
    metadata = render_alignment_svg({"ref": "ACGT", "sample": "TGCA"}, "ref", output)
    marks = events(output, "sample")
    assert metadata["comparison_totals"]["snp"] == 4
    assert {mark.tag for mark in marks} == {f"{SVG}line"}
    assert {mark.get("data-event") for mark in marks} == {"snp"}
    assert len({mark.get("stroke") for mark in marks}) == 4
    assert [mark.get("data-base") for mark in marks] == list("TGCA")


def test_full_alignment_preserves_each_insertion_column_and_angle(tmp_path: Path) -> None:
    output = tmp_path / "indels.svg"
    metadata = render_alignment_svg({"ref": "A-CG-T", "sample": "AT-GCT"}, "ref", output)
    marks = {int(mark.get("data-column")): mark for mark in events(output, "sample")}
    assert {column: mark.get("data-event") for column, mark in marks.items()} == {
        1: "insertion",
        2: "deletion",
        4: "insertion",
    }
    assert marks[1].get("data-angle-deg") == "60.00000000"
    assert marks[4].get("data-angle-deg") == "240.00000000"
    assert marks[1].get("data-reference-position") == "1"
    assert marks[4].get("data-reference-position") == "3"
    assert marks[1].tag == f"{SVG}polygon"
    assert marks[2].tag == f"{SVG}path"
    assert metadata["alignment_columns"] == 6
    assert metadata["plotted_columns"] == list(range(6))


def test_excluding_full_column_leaves_angle_slot_blank(tmp_path: Path) -> None:
    output = tmp_path / "excluded.svg"
    metadata = render_alignment_svg(
        {"ref": "A-CG-T", "sample": "AT-GCT"},
        "ref",
        output,
        excluded_columns={1},
    )
    marks = {int(mark.get("data-column")): mark for mark in events(output, "sample")}
    assert set(marks) == {2, 4}
    assert marks[4].get("data-angle-deg") == "240.00000000"
    assert metadata["angular_slots"] == 6
    assert metadata["columns"] == 5


def test_variable_mode_uses_selected_rows_and_preserves_original_column_ids(tmp_path: Path) -> None:
    records = {"ref": "A-CN", "first": "ATCN", "unselected": "AGTN"}
    output = tmp_path / "selected.svg"
    selected = render_alignment_svg(
        records,
        "ref",
        output,
        selected_ids=["first"],
        coordinate_mode="variable",
    )
    assert selected["selected_ids"] == ["ref", "first"]
    assert selected["plotted_columns"] == [1]
    assert selected["columns"] == selected["angular_slots"] == 1
    assert events(output, "first")[0].get("data-column") == "1"
    assert events(output, "first")[0].get("data-angle-deg") == "0.00000000"
    root = ET.parse(output).getroot()
    assert "Compressed variable mode" in root.find(f"{SVG}desc").text
    all_rows = render_alignment_svg(
        records,
        "ref",
        tmp_path / "all.svg",
        coordinate_mode="variable",
    )
    assert all_rows["plotted_columns"] == [1, 2]


def test_known_sample_difference_can_be_informative_with_unknown_reference(tmp_path: Path) -> None:
    output = tmp_path / "unknown-ref-variable.svg"
    metadata = render_alignment_svg(
        {"ref": "N", "first": "A", "second": "C"},
        "ref",
        output,
        coordinate_mode="variable",
    )
    assert metadata["plotted_columns"] == [0]
    assert metadata["comparison_totals"]["snp"] == 0
    assert {mark.get("data-event") for mark in events(output, "first")} == {"unknown"}


def test_ids_are_xml_text_and_reference_is_always_inner_ring(tmp_path: Path) -> None:
    malicious = '<script>alert("x")</script>&accession'
    reference = 'ref"&<one>'
    output = tmp_path / "escaped.svg"
    metadata = render_alignment_svg(
        {malicious: "AT", "other": "AG", reference: "AC"},
        reference,
        output,
        selected_ids=["other", reference, malicious],
        diameter_mm=320,
    )
    root = ET.parse(output).getroot()
    assert root.get("width") == "320mm"
    assert not list(root.iter(f"{SVG}script"))
    groups = [group for group in root.iter(f"{SVG}g") if group.get("data-record-id") is not None]
    assert [group.get("data-record-id") for group in groups] == [reference, "other", malicious]
    assert float(groups[0].get("data-radius")) < float(groups[1].get("data-radius"))
    embedded = json.loads(root.find(f"{SVG}metadata").text)
    assert embedded == metadata
    assert metadata["selected_ids"] == [reference, "other", malicious]


@pytest.mark.parametrize(
    ("records", "reference", "options", "message"),
    [
        ({"ref": "AC", "short": "A"}, "ref", {}, "equal lengths"),
        ({"ref": "AC"}, "absent", {}, "Unknown reference ID"),
        ({"ref": "AC"}, "ref", {"selected_ids": ["absent"]}, "Unknown selected"),
        ({"ref": "AC"}, "ref", {"selected_ids": ["ref", "ref"]}, "duplicate"),
        ({"ref": "AC"}, "ref", {"excluded_columns": [2]}, "outside the alignment"),
        ({"ref": "AC"}, "ref", {"excluded_columns": [True]}, "integer alignment"),
        ({"ref": "AC"}, "ref", {"excluded_columns": [0, 1]}, "No alignment columns"),
        ({"ref": "AX"}, "ref", {}, "invalid DNA symbols"),
        ({"ref": "--"}, "ref", {}, "entirely gaps"),
        ({"ref": "AC"}, "ref", {"diameter_mm": float("nan")}, "finite positive"),
    ],
)
def test_invalid_inputs_raise_informative_errors(
    tmp_path: Path,
    records: dict[str, str],
    reference: str,
    options: dict[str, object],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        render_alignment_svg(records, reference, tmp_path / "invalid.svg", **options)
