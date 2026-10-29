"""Draw a geographical provenance-root concept without inventing sample places.

Run ``uv run --group preview python scripts/build_root_map.py``. Archaeological
endpoints come directly from data/sites.csv. The exhibition symbol is a
provisional diagrammatic placement in the Matsu area, not an approved venue.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from xml.etree import ElementTree as ET

import build_sample_map as base

ROOT = Path(__file__).resolve().parents[1]
SVG_NS = base.SVG_NS
BROWN = "#795033"
LIGHT_BROWN = "#a77850"
STUMP_LONGITUDE, STUMP_LATITUDE = 119.95, 26.05


def roots_path(x: float, y: float, site: str) -> str:
    start_x, start_y = base.project(STUMP_LONGITUDE, STUMP_LATITUDE)
    controls = {
        "Liangdao": (start_x + 15, start_y - 95, x + 33, y - 66),
        "Suogang": (start_x - 100, start_y + 140, x - 110, y - 120),
        "Hanben": (start_x + 105, start_y + 115, x - 105, y - 168),
        "Gongguan": (start_x + 240, start_y + 25, x + 270, y - 210),
    }
    first_x, first_y, second_x, second_y = controls[site]
    return (
        f"M {start_x:.8f},{start_y:.8f} "
        f"C {first_x:.8f},{first_y:.8f} {second_x:.8f},{second_y:.8f} {x:.8f},{y:.8f}"
    )


def draw_roots(parent: ET.Element, sites: list[dict[str, str]]) -> None:
    for row in sorted(sites, key=lambda item: int(item["ancient_individuals"]), reverse=True):
        x, y = base.project(float(row["longitude"]), float(row["latitude"]))
        count = int(row["ancient_individuals"])
        width = 1.8 + 0.9 * math.sqrt(count)
        curve = base.sub(
            parent,
            "path",
            d=roots_path(x, y, row["site"]),
            fill="none",
            stroke=BROWN,
            **{
                "stroke-width": f"{width:.4f}",
                "stroke-linecap": "round",
                "data-site": row["site"],
                "data-latitude": row["latitude"],
                "data-longitude": row["longitude"],
                "data-count": str(count),
                "data-endpoint-x": f"{x:.8f}",
                "data-endpoint-y": f"{y:.8f}",
            },
        )
        base.sub(curve, "title").text = (
            f"Designed provenance root to {row['site']}: {row['latitude']} N, "
            f"{row['longitude']} E. Endpoint is the CSV-projected archaeological site. "
            "The curve is not a migration route, measured root, or phylogenetic branch."
        )


def cap(
    parent: ET.Element,
    x: float,
    y: float,
    radius: float,
    y_radius: float,
    holes: list[dict[str, str]],
) -> None:
    base.sub(
        parent,
        "ellipse",
        cx=str(x),
        cy=str(y),
        rx=str(radius),
        ry=str(y_radius),
        fill="#d7dddf",
        stroke="#64767e",
        **{"stroke-width": "1.2"},
    )
    if not holes:
        return
    for radius_mm in sorted({float(row["radius_mm"]) for row in holes}):
        fraction = radius_mm / 250.0
        base.sub(
            parent,
            "ellipse",
            cx=str(x),
            cy=str(y),
            rx=f"{radius * fraction:.4f}",
            ry=f"{y_radius * fraction:.4f}",
            fill="none",
            stroke="#89989e",
            **{"stroke-width": "0.35"},
        )
    for row in holes:
        hole_x = x + (float(row["x_mm"]) - 250.0) / 250.0 * radius
        hole_y = y + (float(row["y_mm"]) - 250.0) / 250.0 * y_radius
        diameter = float(row["hole_diameter_mm"])
        base.sub(
            parent,
            "ellipse",
            cx=f"{hole_x:.4f}",
            cy=f"{hole_y:.4f}",
            rx=f"{diameter / 500 * radius:.4f}",
            ry=f"{diameter / 500 * y_radius:.4f}",
            fill="#25343a",
        )
    base.sub(
        parent,
        "ellipse",
        cx=str(x),
        cy=str(y),
        rx=f"{radius * 30 / 250:.4f}",
        ry=f"{y_radius * 30 / 250:.4f}",
        fill="none",
        stroke="#89989e",
        **{"stroke-width": "0.4"},
    )


def stump_symbol(parent: ET.Element, holes: list[dict[str, str]]) -> None:
    x, y = base.project(STUMP_LONGITUDE, STUMP_LATITUDE)
    symbol = base.sub(parent, "g", id="provisional-exhibition-stump")
    base.sub(
        symbol, "title"
    ).text = "Diagrammatic exhibition stump in the Matsu area. Exact venue unknown; symbol not to map scale."
    base.sub(
        symbol,
        "circle",
        cx=str(x),
        cy=str(y),
        r="46",
        fill="#af8356",
        stroke=BROWN,
        **{"stroke-width": "2"},
    )
    for angle in range(0, 360, 28):
        radians = math.radians(angle)
        x1, y1 = x + 40 * math.cos(radians), y + 40 * math.sin(radians)
        x2, y2 = x + 45 * math.cos(radians + 0.04), y + 45 * math.sin(radians + 0.04)
        base.sub(
            symbol,
            "line",
            x1=str(x1),
            y1=str(y1),
            x2=str(x2),
            y2=str(y2),
            stroke="#705037",
            **{"stroke-width": "1.5"},
        )
    cap(symbol, x, y, 38, 38, holes)
    base.sub(
        parent,
        "path",
        d=f"M {x - 42:.3f},{y + 18:.3f} L 363,392",
        fill="none",
        stroke="#6c7778",
        **{"stroke-width": "1", "stroke-dasharray": "4 4"},
    )
    base.sub(
        parent,
        "rect",
        x="112",
        y="357",
        width="255",
        height="101",
        rx="7",
        fill="#ffffff",
        stroke="#d1dcda",
    )
    base.text(parent, "Stump of Life", 127, 385, 25, BROWN, "bold")
    base.text(parent, "Matsu exhibition area · provisional", 127, 412, 15)
    base.text(parent, "Exact venue to be agreed", 127, 437, 14, base.MUTED)


def endpoint_label(parent: ET.Element, row: dict[str, str]) -> None:
    x, y = base.project(float(row["longitude"]), float(row["latitude"]))
    name = row["site"]
    positions = {
        "Liangdao": (610.0, 207.0, 355.0, "Matsu archipelago"),
        "Suogang": (120.0, 817.0, 292.0, "Penghu archipelago"),
        "Hanben": (810.0, 430.0, 245.0, "Blihun site, Taiwan"),
        "Gongguan": (790.0, 996.0, 265.0, "Green Island"),
    }
    left, top, width, region = positions[name]
    edge_x = left if x < left else left + width
    edge_y = min(top + 122, max(top + 22, y))
    base.sub(
        parent,
        "line",
        x1=str(x),
        y1=str(y),
        x2=str(edge_x),
        y2=str(edge_y),
        stroke=base.ACCENT,
        **{"stroke-width": "1.3"},
    )
    marker = base.sub(
        parent,
        "circle",
        cx=f"{x:.8f}",
        cy=f"{y:.8f}",
        r="6.5",
        fill=base.ACCENT,
        stroke="#ffffff",
        **{
            "stroke-width": "2",
            "data-site-endpoint": name,
            "data-latitude": row["latitude"],
            "data-longitude": row["longitude"],
        },
    )
    base.sub(marker, "title").text = f"{name}; {row['precision_note']}; {row['coordinate_source']}"
    base.sub(
        parent,
        "rect",
        x=str(left),
        y=str(top),
        width=str(width),
        height="140",
        rx="7",
        fill="#ffffff",
        stroke="#d1dcda",
    )
    base.text(
        parent,
        f"{name} · {row['ancient_individuals']}",
        left + 15,
        top + 29,
        25,
        base.ACCENT,
        "bold",
    )
    base.text(parent, region, left + 15, top + 54, 17)
    base.text(parent, f"{row['latitude']}° N", left + 15, top + 79, 16)
    base.text(parent, f"{row['longitude']}° E", left + 15, top + 100, 16)
    precision = (
        "Rounded published coordinate" if name == "Suogang" else "Archaeological-site coordinate"
    )
    base.text(parent, precision, left + 15, top + 124, 12.5, base.MUTED)


def elevation(
    parent: ET.Element,
    geometry: dict[str, object],
    sites: list[dict[str, str]],
    holes: list[dict[str, str]],
) -> None:
    base.text(parent, "WOOD STUMP + METAL DNA CAP", 1120, 804, 17, BROWN, "bold")
    base.text(parent, "Material / elevation concept", 1120, 835, 18, base.MUTED)
    # A perspective regional-map relief carries the same verified coast and roots.
    a, b = 310 / base.MAP_W, -70 / base.MAP_W
    c, d = 90 / base.MAP_H, 85 / base.MAP_H
    e, f = 1125 - a * base.MAP_X - c * base.MAP_Y, 1035 - b * base.MAP_X - d * base.MAP_Y
    base.sub(
        parent,
        "polygon",
        points="1125,1035 1435,965 1525,1050 1215,1120 1215,1130 1525,1060 1525,1050 1435,975 1125,1045",
        fill="#b9b09f",
        stroke="#9f9583",
    )
    mini = base.sub(
        parent,
        "g",
        transform=f"matrix({a:.8f} {b:.8f} {c:.8f} {d:.8f} {e:.8f} {f:.8f})",
        **{"clip-path": "url(#root-map-clip)"},
    )
    base.sub(
        mini,
        "rect",
        x=str(base.MAP_X),
        y=str(base.MAP_Y),
        width=str(base.MAP_W),
        height=str(base.MAP_H),
        fill=base.OCEAN,
    )
    base.draw_land(mini, geometry)
    draw_roots(mini, sites)
    map_x, map_y = base.project(STUMP_LONGITUDE, STUMP_LATITUDE)
    x, floor_y = a * map_x + c * map_y + e, b * map_x + d * map_y + f
    top = floor_y - 143
    body = (
        f"M {x - 49:.2f},{top:.2f} Q {x - 57:.2f},{top + 35:.2f} {x - 52:.2f},{top + 75:.2f} "
        f"L {x - 57:.2f},{floor_y - 9:.2f} Q {x:.2f},{floor_y + 21:.2f} {x + 57:.2f},{floor_y - 8:.2f} "
        f"L {x + 51:.2f},{top + 76:.2f} Q {x + 56:.2f},{top + 35:.2f} {x + 50:.2f},{top:.2f} Z"
    )
    base.sub(parent, "path", d=body, fill="#9a7049", stroke=BROWN, **{"stroke-width": "1.5"})
    for index in range(9):
        offset = -43 + index * 10.5
        base.sub(
            parent,
            "path",
            d=f"M {x + offset:.2f},{top + 9:.2f} Q {x + offset - 4:.2f},{top + 64:.2f} {x + offset + 2:.2f},{floor_y - 3:.2f}",
            fill="none",
            stroke="#775338" if index % 2 else "#c09665",
            **{"stroke-width": "2.5" if index % 2 else "1.5"},
        )
    cap(parent, x, top, 62, 19, holes)
    base.sub(
        parent,
        "path",
        d=f"M {x + 61:.2f},{top:.2f} L 1420,861",
        fill="none",
        stroke=base.MUTED,
        **{"stroke-width": "1"},
    )
    base.text(parent, "Metal cap", 1422, 859, 15)
    base.sub(
        parent,
        "path",
        d=f"M {x + 52:.2f},{top + 83:.2f} L 1420,947",
        fill="none",
        stroke=base.MUTED,
        **{"stroke-width": "1"},
    )
    base.text(parent, "Natural wood", 1422, 946, 14)
    base.text(parent, "Regional map relief", 1245, 1172, 17, base.MUTED)
    base.text(parent, "Schematic materials; symbol not to scale", 1120, 1213, 15, base.MUTED)


def main() -> None:
    with (ROOT / "data" / "sites.csv").open(newline="") as handle:
        sites = list(csv.DictReader(handle))
    geometry, geography_provenance = base.land_geometry()
    holes_path = ROOT / "output" / "taiwan" / "metal" / "holes.csv"
    holes: list[dict[str, str]] = []
    if holes_path.exists():
        with holes_path.open(newline="") as handle:
            holes = list(csv.DictReader(handle))
    root = ET.Element(
        f"{{{SVG_NS}}}svg",
        {
            "width": "1550",
            "height": "1420",
            "viewBox": "0 0 1550 1420",
            "role": "img",
            "aria-labelledby": "root-title root-description",
            "font-family": "DejaVu Sans, Arial, sans-serif",
        },
    )
    base.sub(
        root, "title", id="root-title"
    ).text = "Stump of Life: roots to the archaeological source sites"
    base.sub(root, "desc", id="root-description").text = (
        "Exhibition concept: a provisional stump symbol in the Matsu area connects by designed "
        "provenance roots to the exact CSV-projected archaeological-site coordinates for Liangdao, "
        "Suogang, Hanben and Gongguan. One Hanben endpoint represents 45 individuals sharing a site "
        "coordinate, and Gongguan represents one. Curves depict provenance, not migration or phylogeny. "
        "The 550 modern Taiwan samples have no verified collection coordinates here and are not mapped."
    )
    base.sub(root, "metadata").text = json.dumps(
        {
            "sites": sites,
            "base_map": geography_provenance,
            "exhibition_position": {
                "longitude": STUMP_LONGITUDE,
                "latitude": STUMP_LATITUDE,
                "meaning": "diagrammatic provisional Matsu-area symbol; no venue claimed",
            },
            "root_width_formula": "1.8 + 0.9 * sqrt(site individual count), SVG pixels",
            "root_geometry": "designed cubic curves; endpoints equal projected data/sites.csv coordinates",
            "cap_hole_source": str(holes_path.relative_to(ROOT)) if holes else None,
            "cap_hole_count": len(holes),
            "modern_collection_points_plotted": False,
        },
        ensure_ascii=False,
    )
    base.sub(root, "rect", width="1550", height="1420", fill="#fcfbf7")
    base.text(root, "Stump of Life · roots to the source sites", 78, 61, 36, weight="bold")
    base.text(
        root, "A geographical provenance concept for an exhibition in Matsu", 80, 99, 22, base.MUTED
    )
    base.text(
        root,
        "Four archaeological endpoints · 50 ancient individuals · one proposed exhibition stump",
        80,
        133,
        19,
    )
    defs = base.sub(root, "defs")
    clip = base.sub(defs, "clipPath", id="root-map-clip", clipPathUnits="userSpaceOnUse")
    base.sub(
        clip,
        "rect",
        x=str(base.MAP_X),
        y=str(base.MAP_Y),
        width=str(base.MAP_W),
        height=str(base.MAP_H),
    )
    map_group = base.sub(root, "g", **{"clip-path": "url(#root-map-clip)"})
    base.sub(
        map_group,
        "rect",
        x=str(base.MAP_X),
        y=str(base.MAP_Y),
        width=str(base.MAP_W),
        height=str(base.MAP_H),
        fill=base.OCEAN,
    )
    for longitude in range(118, 124):
        x, _ = base.project(longitude, base.NORTH)
        base.sub(
            map_group,
            "line",
            x1=str(x),
            y1=str(base.MAP_Y),
            x2=str(x),
            y2=str(base.MAP_Y + base.MAP_H),
            stroke="#d1e0e1",
            **{"stroke-width": "0.8"},
        )
        base.text(
            root,
            f"{longitude}° E",
            x,
            base.MAP_Y + base.MAP_H + 28,
            15,
            base.MUTED,
            anchor="middle",
        )
    for latitude in range(21, 28):
        _, y = base.project(base.WEST, latitude)
        base.sub(
            map_group,
            "line",
            x1=str(base.MAP_X),
            y1=str(y),
            x2=str(base.MAP_X + base.MAP_W),
            y2=str(y),
            stroke="#d1e0e1",
            **{"stroke-width": "0.8"},
        )
        base.text(root, f"{latitude}° N", base.MAP_X - 12, y + 5, 15, base.MUTED, anchor="end")
    base.draw_land(map_group, geometry)
    base.text(
        map_group, "Fujian coast", *base.project(118.35, 25.75), 23, base.MUTED, anchor="middle"
    )
    base.text(map_group, "Taiwan", *base.project(120.95, 23.85), 30, base.MUTED, "bold", "middle")
    base.text(
        map_group, "Taiwan Strait", *base.project(119.30, 24.45), 20, "#79979c", anchor="middle"
    )
    draw_roots(map_group, sites)
    stump_symbol(root, holes)
    for row in sites:
        endpoint_label(root, row)
    base.sub(
        root,
        "rect",
        x=str(base.MAP_X),
        y=str(base.MAP_Y),
        width=str(base.MAP_W),
        height=str(base.MAP_H),
        fill="none",
        stroke="#adc3c5",
    )

    x = 1120
    base.text(root, "A GEOGRAPHICAL SCULPTURE", x, 207, 17, BROWN, "bold")
    for index, line in enumerate(
        [
            "One stump in the Matsu exhibition area.",
            "Roots end at the source sites shown",
            "in the archaeological metadata.",
            "Exact venue remains to be agreed.",
        ]
    ):
        base.text(root, line, x, 245 + index * 28, 17)
    for index, line in enumerate(
        [
            "Thicker roots indicate more individuals.",
            "45 people at Hanben share one endpoint;",
            "Gongguan has one individual.",
        ]
    ):
        base.text(root, line, x, 385 + index * 27, 16, base.MUTED)
    base.sub(root, "line", x1=str(x), y1="487", x2="1495", y2="487", stroke="#d2dcda")
    base.text(root, "MODERN COMPARISON PANEL", x, 525, 17, base.GREEN, "bold")
    base.text(root, "550", x, 583, 51, base.GREEN, "bold")
    base.text(root, "Taiwan participants · 12 study groups", x, 615, 18)
    base.text(root, "Ko et al. (2014)", x, 645, 17, base.MUTED)
    for index, line in enumerate(
        [
            "Collection coordinates are unverified here;",
            "no modern collection-point roots are drawn.",
        ]
    ):
        base.text(root, line, x, 687 + index * 25, 15.5, base.MUTED)
    elevation(root, geometry, sites, holes)
    base.text(
        root,
        "Designed roots show geographical provenance; they do not represent migration routes, genealogical branches or measured travel paths.",
        80,
        1321,
        16,
        base.MUTED,
    )
    base.text(
        root,
        "Archaeological-site coordinates are not individual burial GPS. Suogang remains at its published rounded coordinate; Liangdao uses the official heritage location.",
        80,
        1351,
        14.5,
        base.MUTED,
    )
    base.text(
        root,
        "Sources: data/sites.csv · Ministry of Culture; Yang et al. 2020; Wang et al. 2021 · Natural Earth 1:50m public-domain coast.",
        80,
        1380,
        15,
        base.MUTED,
    )
    output = ROOT / "docs" / "images"
    output.mkdir(parents=True, exist_ok=True)
    svg_path = output / "root-map.svg"
    png_path = output / "root-map.png"
    ET.indent(root, space="  ")
    ET.ElementTree(root).write(svg_path, encoding="utf-8", xml_declaration=True)
    import cairosvg

    cairosvg.svg2png(url=str(svg_path), write_to=str(png_path), output_width=1550)
    print(
        json.dumps(
            {
                "svg": str(svg_path),
                "png": str(png_path),
                "cap_holes": len(holes),
                "root_endpoints": len(sites),
            }
        )
    )


if __name__ == "__main__":
    main()
