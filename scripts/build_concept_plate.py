"""Draw the dated sample plan from actual aligned variants, without invented ancient DNA.

Run after the example alignment and metal export. Output is an annotated design
plan: reserved ancient tracks carry no fabricated mismatch pattern.
"""

import csv
import json
from pathlib import Path
from xml.etree import ElementTree as ET

import typer

from stump_of_life.io import read_fasta
from stump_of_life.metal import export_perforation

app = typer.Typer()
SVG = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG)


def element(parent: ET.Element, tag: str, **attrs: object) -> ET.Element:
    return ET.SubElement(parent, f"{{{SVG}}}{tag}", {key: str(value) for key, value in attrs.items()})


def text(parent: ET.Element, value: str, x: float, y: float, size: float = 20,
         fill: str = "#243a3b", **attrs: object) -> None:
    node = element(parent, "text", x=x, y=y, fill=fill, **{"font-size": size, **attrs})
    node.text = value


@app.command()
def build(
    alignment: Path = Path("output/taiwan/aligned.fasta"),
    output: Path = Path("docs/images/dated-sample-plan.svg"),
) -> None:
    """Generate the 500 mm cap, dated tracks and a real-trunk mounting section."""
    records = read_fasta(alignment, aligned=True)
    metal_output = Path("output/concept/metal")
    geometry = export_perforation(records, "KF540505.1", metal_output)
    with (metal_output / "holes.csv").open(newline="") as handle:
        holes = list(csv.DictReader(handle))
    root = ET.Element(f"{{{SVG}}}svg", {
        "width": "1800", "height": "1320", "viewBox": "0 0 1800 1320",
        "font-family": "DejaVu Sans, Arial, sans-serif",
    })
    element(root, "rect", x=0, y=0, width=1800, height=1320, fill="#f5f2e9")
    text(root, "STUMP OF LIFE", 70, 78, 44, **{"font-weight": "bold"})
    text(root, "Real tree trunk + replaceable mitochondrial DNA cross-section", 72, 120, 25)
    text(root, "Dated sample placement / 500 mm top / 5 mm substitution holes", 72, 156, 19, "#63726b")

    cx, cy, scale = 495.0, 650.0, 1.62
    radius = 250 * scale
    element(root, "circle", cx=cx, cy=cy, r=radius + 19, fill="#826747")
    element(root, "circle", cx=cx, cy=cy, r=radius, fill="#d8ddd8", stroke="#4d625c", **{"stroke-width": 2})
    # These guides have declared proposed radii; missing sequences receive no holes.
    pending = [
        (55, "Liangdao 2", "#4b8279"),
        (92, "Suogang B1", "#b18b50"), (100, "Suogang B3", "#b18b50"),
        (125, "Gongguan I13721", "#a9755b"),
        (145, "Hanben I3620", "#7a8064"), (160, "Hanben I8081", "#7a8064"),
    ]
    for row_radius, name, color in pending:
        ring = element(root, "circle", cx=cx, cy=cy, r=row_radius * scale,
                       fill="none", stroke=color, **{"stroke-width": 1.6, "stroke-dasharray": "6 6"})
        element(ring, "title").text = f"{name}: reserved track; no DNA differences invented"
    element(root, "circle", cx=cx, cy=cy, r=30 * scale, fill="#c6d8cd", stroke="#326d62", **{"stroke-width": 2})
    for guide in geometry["guides"][1:]:
        element(root, "circle", cx=cx, cy=cy, r=guide["radius_mm"] * scale,
                fill="none", stroke="#82968d", **{"stroke-width": 0.8})
    for hole in holes:
        mark = element(root, "circle", cx=cx + (float(hole["x_mm"]) - 250) * scale,
                       cy=cy + (float(hole["y_mm"]) - 250) * scale, r=2.5 * scale, fill="#253d39")
        element(mark, "title").text = (
            f"{hole['sample_id']}: reference position {hole['reference_position']}; "
            f"{hole['reference_base']} to {hole['sample_base']}"
        )
    text(root, "Liangdao 1", cx, cy - 8, 14, **{"text-anchor": "middle", "font-weight": "bold"})
    text(root, "reference", cx, cy + 14, 12, **{"text-anchor": "middle"})
    # Ring numbers sit in a clear guide ray. They refer to the timeline panel.
    for number, row_radius, color in [("2", 55, "#4b8279"), ("3", 96, "#b18b50"),
                                      ("4", 125, "#a9755b"), ("5", 152.5, "#7a8064"),
                                      ("6", 209, "#476b78")]:
        x, y = cx - row_radius * scale, cy
        element(root, "circle", cx=x, cy=y, r=15, fill="#f5f2e9", stroke=color)
        text(root, number, x, y + 5, 15, color, **{"text-anchor": "middle", "font-weight": "bold"})
    # Physical dimension line, independent of the artistic period-band scale.
    for x in (cx - radius, cx + radius):
        element(root, "line", x1=x, y1=cy + radius + 30, x2=x,
                y2=cy + radius + 55, stroke="#52685f", **{"stroke-width": 1.2})
    element(root, "line", x1=cx - radius, y1=cy + radius + 43,
            x2=cx + radius, y2=cy + radius + 43, stroke="#52685f")
    text(root, "500 mm DNA insert", cx, cy + radius + 76, 21, **{"text-anchor": "middle"})
    text(root, "Chronology increases outwards; radii use period bands, not linear time.", 72, 1180, 19)
    text(root, "Within a shared date band, track order does not represent generations.", 72, 1210, 19)
    text(root, "Solid outer tracks and holes: real sequence data. Dashed tracks: proposed positions.", 72, 1240, 19)

    timeline_x = 1000
    text(root, "CHRONOLOGY AND SAMPLE POSITIONS", timeline_x, 220, 24, **{"font-weight": "bold"})
    timeline = [
        ("1", "Liangdao 1 / Matsu", "8,060–8,320 cal BP", "r = 30 mm · KF540505.1 · complete reference", "#326d62"),
        ("2", "Liangdao 2 / Matsu", "7,510–7,590 cal BP", "r = 55 mm · sequence validation pending", "#4b8279"),
        ("3", "Suogang B1 and B3 / Penghu", "4,793–4,287 cal BP · associated shells", "r = 92, 100 mm · reserved positions", "#b18b50"),
        ("4", "Gongguan / Green Island", "I13721 · 1366–1126 cal BCE", "r = 125 mm · reserved position", "#a9755b"),
        ("5", "Hanben / Taiwan", "I3620: 22–201 CE · I8081: 439–586 CE", "r = 145, 160 mm · 2 selected from 45 individuals", "#7a8064"),
        ("6", "Modern Taiwanese sequences", "Collected 1998–2001", "9 selected from 550 · r = 175–240 mm", "#476b78"),
    ]
    for index, (number, name, date, placement, color) in enumerate(timeline):
        y = 280 + index * 110
        element(root, "circle", cx=timeline_x + 14, cy=y - 5, r=16, fill=color)
        text(root, number, timeline_x + 14, y + 1, 16, "#ffffff", **{"text-anchor": "middle"})
        text(root, name, timeline_x + 46, y, 22, **{"font-weight": "bold"})
        text(root, date, timeline_x + 46, y + 29, 18)
        text(root, placement, timeline_x + 46, y + 55, 17, "#6c7569")

    # An honest section drawing, not a fabricated photorealistic material claim.
    text(root, "FITTED CAP ON A REAL LOG", timeline_x, 975, 23, **{"font-weight": "bold"})
    element(root, "rect", x=1035, y=1050, width=355, height=170, rx=10,
            fill="#a3835b", stroke="#6a5539", **{"stroke-width": 2})
    for x in (1060, 1110, 1155, 1218, 1260, 1325, 1370):
        element(root, "path", d=f"M {x} 1060 Q {x - 12} 1120 {x + 4} 1215",
                stroke="#806243", fill="none", **{"stroke-width": 2})
    element(root, "rect", x=1026, y=1027, width=373, height=9, rx=3, fill="#a7b5b0", stroke="#49635d")
    for x in (1060, 1364):
        element(root, "rect", x=x, y=1036, width=14, height=14, fill="#4d625c")
    text(root, "Replaceable metal top", 1430, 1040, 17)
    text(root, "Ventilated gap + supports", 1430, 1070, 17)
    text(root, "Seasoned log body", 1430, 1110, 17)
    text(root, "Target height: 43–48 cm", 1430, 1140, 17)
    text(root, "Mounting detail requires", 1430, 1190, 16, "#6c7569")
    text(root, "fabricator confirmation.", 1430, 1214, 16, "#6c7569")

    footer = (
        f"Design plan: 16 tracks, including 6 reserved ancient tracks; {geometry['hole_count']} real modern substitution holes. "
        "Not a structural or fabrication approval."
    )
    text(root, footer, 72, 1282, 15, "#6b7167")
    element(root, "metadata").text = json.dumps({
        "reference": "KF540505.1", "real_holes": geometry["hole_count"],
        "reserved_ancient_tracks": [name for _, name, _ in pending],
        "all_ancient_individuals_in_scope": 50,
        "all_modern_individuals_in_scope": 550,
        "sources": ["https://doi.org/10.1073/pnas.2026132118",
                    "https://doi.org/10.1126/science.aba0909",
                    "https://doi.org/10.1038/s41586-021-03336-2",
                    "https://doi.org/10.1016/j.ajhg.2014.02.003"],
        "suogang_dates": "Associated shells, not direct skeletal dates; within-band order is arbitrary.",
        "geometry": geometry,
    }, ensure_ascii=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(root)
    ET.ElementTree(root).write(output, encoding="utf-8", xml_declaration=True)
    typer.echo(str(output.resolve()))


if __name__ == "__main__":
    app()
