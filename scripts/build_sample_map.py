"""Build an offline geographical SVG and PNG from verified sample-site metadata.

Run ``uv run --group preview python scripts/build_sample_map.py``. The first run
downloads public-domain Natural Earth 1:50m land geometry into data/geography;
later runs use that cached file. Modern population names are a study legend,
not invented individual collection coordinates or migration routes.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
GEOGRAPHY = ROOT / "data" / "geography"
OUTPUT = ROOT / "docs" / "images"
LAND_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
    "geojson/ne_50m_land.geojson"
)
SVG_NS = "http://www.w3.org/2000/svg"
WEST, EAST, SOUTH, NORTH = 117.5, 123.5, 21.0, 27.0
MAP_X, MAP_Y, MAP_W = 80.0, 170.0, 990.0
MAP_H = MAP_W / math.cos(math.radians(24.0))
WIDTH, HEIGHT = 1550, 1420
INK = "#1e3239"
MUTED = "#60747a"
OCEAN = "#e8f1f2"
LAND = "#e9e7dc"
ACCENT = "#b9472c"
GREEN = "#136b61"
GROUPS = [
    "Amis",
    "Atayal",
    "Bunun",
    "Hakka",
    "Minnan",
    "Paiwan",
    "Makatao",
    "Puyuma",
    "Rukai",
    "Saisiyat",
    "Tao",
    "Tsou",
]
ET.register_namespace("", SVG_NS)


def sub(parent: ET.Element, tag: str, **attributes: str) -> ET.Element:
    return ET.SubElement(parent, f"{{{SVG_NS}}}{tag}", attributes)


def text(
    parent: ET.Element,
    value: str,
    x: float,
    y: float,
    size: float = 20.0,
    color: str = INK,
    weight: str = "normal",
    anchor: str = "start",
) -> ET.Element:
    item = sub(
        parent,
        "text",
        x=f"{x:.2f}",
        y=f"{y:.2f}",
        fill=color,
        **{"font-size": f"{size:g}", "font-weight": weight, "text-anchor": anchor},
    )
    item.text = value
    return item


def project(longitude: float, latitude: float) -> tuple[float, float]:
    return (
        MAP_X + (longitude - WEST) / (EAST - WEST) * MAP_W,
        MAP_Y + (NORTH - latitude) / (NORTH - SOUTH) * MAP_H,
    )


def land_geometry() -> tuple[dict[str, Any], dict[str, object]]:
    GEOGRAPHY.mkdir(parents=True, exist_ok=True)
    path = GEOGRAPHY / "ne_50m_land.geojson"
    provenance_path = GEOGRAPHY / "natural_earth_provenance.json"
    if not path.exists():
        request = Request(LAND_URL, headers={"User-Agent": "stump-of-life public map preview"})
        with urlopen(request, timeout=60) as response:
            raw = response.read()
        geometry = json.loads(raw)
        path.write_bytes(raw)
        provenance = {
            "source": LAND_URL,
            "dataset": "Natural Earth 1:50m land",
            "license": "Public domain",
            "license_url": "https://www.naturalearthdata.com/about/terms-of-use/",
            "retrieved_utc": datetime.now(UTC).isoformat(timespec="seconds"),
            "sha256": hashlib.sha256(raw).hexdigest(),
        }
        provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    else:
        raw = path.read_bytes()
        geometry = json.loads(raw)
        provenance = json.loads(provenance_path.read_text())
        if hashlib.sha256(raw).hexdigest() != provenance["sha256"]:
            raise ValueError("Cached Natural Earth geometry differs from its recorded checksum.")
    return geometry, provenance


def draw_land(parent: ET.Element, geometry: dict[str, Any]) -> None:
    for feature in geometry["features"]:
        geom = feature["geometry"]
        polygons = [geom["coordinates"]] if geom["type"] == "Polygon" else geom["coordinates"]
        for polygon in polygons:
            exterior = polygon[0]
            longitudes = [point[0] for point in exterior]
            latitudes = [point[1] for point in exterior]
            if max(longitudes) < WEST or min(longitudes) > EAST:
                continue
            if max(latitudes) < SOUTH or min(latitudes) > NORTH:
                continue
            paths = []
            for ring in polygon:
                points = [project(float(lon), float(lat)) for lon, lat, *_ in ring]
                paths.append("M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in points) + " Z")
            sub(
                parent,
                "path",
                d=" ".join(paths),
                fill=LAND,
                stroke="#879b9a",
                **{"stroke-width": "1.15", "fill-rule": "evenodd", "stroke-linejoin": "round"},
            )


def site_label(parent: ET.Element, row: dict[str, str]) -> None:
    name = row["site"]
    x, y = project(float(row["longitude"]), float(row["latitude"]))
    settings = {
        "Liangdao": (
            610.0,
            225.0,
            355.0,
            "Matsu archipelago",
            "Official archaeological-site location",
        ),
        "Suogang": (125.0, 815.0, 285.0, "Penghu archipelago", "Rounded study-site coordinate"),
        "Hanben": (
            810.0,
            565.0,
            245.0,
            "Blihun site, Taiwan",
            "Published site coordinate",
        ),
        "Gongguan": (
            790.0,
            1000.0,
            265.0,
            "Green Island",
            "Published site coordinate",
        ),
    }
    left, top, width, region, note = settings[name]
    edge_x = left if x < left else left + width
    edge_y = min(top + 77, max(top + 20, y))
    sub(
        parent,
        "line",
        x1=f"{x:.2f}",
        y1=f"{y:.2f}",
        x2=f"{edge_x:.2f}",
        y2=f"{edge_y:.2f}",
        stroke=ACCENT,
        **{"stroke-width": "1.6"},
    )
    count = int(row["ancient_individuals"])
    radius = 10.0 if count > 10 else 7.0
    marker = sub(
        parent,
        "circle",
        cx=f"{x:.2f}",
        cy=f"{y:.2f}",
        r=f"{radius:g}",
        fill=ACCENT,
        stroke="#ffffff",
        **{"stroke-width": "2.5", "data-site": name},
    )
    sub(marker, "title").text = (
        f"{name}: {row['latitude']} N, {row['longitude']} E; {count} individuals. "
        f"{row['precision_note']}. Source: {row['coordinate_source']}"
    )
    sub(
        parent,
        "rect",
        x=str(left),
        y=str(top),
        width=str(width),
        height="103",
        rx="7",
        fill="#ffffff",
        stroke="#d1dcda",
        **{"stroke-width": "1"},
    )
    text(parent, f"{name} · {count}", left + 15, top + 29, 25, ACCENT, "bold")
    text(parent, region, left + 15, top + 54, 18)
    # The note is factual precision metadata, not a claim of exact burial GPS.
    text(parent, note, left + 15, top + 78, 12.5, MUTED)


def main() -> None:
    with (ROOT / "data" / "sites.csv").open(newline="") as handle:
        sites = list(csv.DictReader(handle))
    geometry, provenance = land_geometry()
    ancient_count = sum(int(site["ancient_individuals"]) for site in sites)
    root = ET.Element(
        f"{{{SVG_NS}}}svg",
        {
            "width": str(WIDTH),
            "height": str(HEIGHT),
            "viewBox": f"0 0 {WIDTH} {HEIGHT}",
            "role": "img",
            "aria-labelledby": "map-title map-description",
            "font-family": "DejaVu Sans, Arial, sans-serif",
        },
    )
    sub(root, "title", id="map-title").text = "Ancient sites and a modern Taiwan study panel"
    sub(root, "desc", id="map-description").text = (
        "Geographical map of Liangdao in the Matsu archipelago (2 ancient individuals), "
        "Suogang in Penghu (2), Hanben in Taiwan (45), and Gongguan on Green Island (1). "
        "Modern Ko et al. 2014 samples comprise 550 people from 12 reported groups in Taiwan; "
        "unverified collection locations are not plotted. Site points do not locate individual burials."
    )
    sub(root, "metadata").text = json.dumps(
        {
            "sites": sites,
            "modern_count": 550,
            "modern_groups": GROUPS,
            "base_map": provenance,
            "extent": {"west": WEST, "east": EAST, "south": SOUTH, "north": NORTH},
            "projection": "Equirectangular; longitude scale corrected at 24 degrees north",
        },
        ensure_ascii=False,
    )
    sub(root, "rect", width=str(WIDTH), height=str(HEIGHT), fill="#fcfbf7")
    text(root, "Ancient sites, modern Taiwan", 78, 61, 39, weight="bold")
    text(
        root,
        "Archaeological locations and the scope of the present-day comparison panel",
        80,
        99,
        22,
        MUTED,
    )
    text(
        root,
        "50 ancient individuals at 4 sites  ·  550 modern individuals across Taiwan",
        80,
        133,
        20,
    )
    defs = sub(root, "defs")
    clip = sub(defs, "clipPath", id="map-clip")
    sub(clip, "rect", x=str(MAP_X), y=str(MAP_Y), width=str(MAP_W), height=f"{MAP_H:.2f}")
    map_group = sub(root, "g", **{"clip-path": "url(#map-clip)"})
    sub(
        map_group,
        "rect",
        x=str(MAP_X),
        y=str(MAP_Y),
        width=str(MAP_W),
        height=f"{MAP_H:.2f}",
        fill=OCEAN,
    )
    for longitude in range(118, 124):
        x, _ = project(longitude, NORTH)
        sub(
            map_group,
            "line",
            x1=str(x),
            y1=str(MAP_Y),
            x2=str(x),
            y2=str(MAP_Y + MAP_H),
            stroke="#d1e0e1",
            **{"stroke-width": "0.8"},
        )
    for latitude in range(21, 28):
        _, y = project(WEST, latitude)
        sub(
            map_group,
            "line",
            x1=str(MAP_X),
            y1=str(y),
            x2=str(MAP_X + MAP_W),
            y2=str(y),
            stroke="#d1e0e1",
            **{"stroke-width": "0.8"},
        )
    draw_land(map_group, geometry)
    for longitude in range(118, 124):
        x, _ = project(longitude, NORTH)
        text(root, f"{longitude}° E", x, MAP_Y + MAP_H + 28, 15, MUTED, anchor="middle")
    for latitude in range(21, 28):
        _, y = project(WEST, latitude)
        text(root, f"{latitude}° N", MAP_X - 12, y + 5, 15, MUTED, anchor="end")
    sub(
        root,
        "rect",
        x=str(MAP_X),
        y=str(MAP_Y),
        width=str(MAP_W),
        height=f"{MAP_H:.2f}",
        fill="none",
        stroke="#adc3c5",
        **{"stroke-width": "1"},
    )
    text(map_group, "Fujian coast", *project(118.35, 26.0), 24, MUTED, anchor="middle")
    text(map_group, "Taiwan", *project(120.95, 23.85), 31, MUTED, "bold", "middle")
    sea_x, sea_y = project(119.6, 24.55)
    sea_label = text(map_group, "Taiwan Strait", sea_x, sea_y, 24, "#7b9a9f", anchor="middle")
    sea_label.set("transform", f"rotate(-49 {sea_x:.2f} {sea_y:.2f})")
    text(map_group, "Pacific Ocean", *project(122.55, 22.15), 22, "#7b9a9f", anchor="middle")
    for row in sites:
        site_label(root, row)
    # Simple north arrow and local 100 km scale; no speculative migration paths.
    text(root, "N", 1015, 210, 18, anchor="middle")
    sub(root, "path", d="M 1015,218 L 1008,245 L 1015,239 L 1022,245 Z", fill=MUTED)
    scale_pixels = MAP_W / 6 * (100 / (111.32 * math.cos(math.radians(24))))
    scale_x, scale_y = MAP_X + 36, MAP_Y + MAP_H - 42
    sub(
        root,
        "path",
        d=f"M {scale_x},{scale_y - 7} V {scale_y} H {scale_x + scale_pixels} V {scale_y - 7}",
        fill="none",
        stroke=MUTED,
        **{"stroke-width": "2"},
    )
    text(
        root,
        "100 km at 24° N",
        scale_x + scale_pixels / 2,
        scale_y + 23,
        14,
        MUTED,
        anchor="middle",
    )

    panel_x = 1120
    text(root, "MODERN COMPARISON", panel_x, 205, 17, GREEN, "bold")
    text(root, "550", panel_x, 269, 57, GREEN, "bold")
    text(root, "individuals · 12 study groups", panel_x, 302, 20)
    text(root, "Ko et al. · 2014", panel_x, 334, 18, MUTED)
    text(root, "Country reported: Taiwan", panel_x, 370, 18)
    for index, group in enumerate(GROUPS):
        text(root, group, panel_x + (index // 6) * 170, 410 + (index % 6) * 31, 20)
    for index, line in enumerate(
        [
            "Collection points are not geolocated.",
            "Group names describe the study panel;",
            "they are not new map coordinates.",
        ]
    ):
        text(root, line, panel_x, 625 + index * 26, 16, MUTED)
    sub(root, "line", x1=str(panel_x), y1="721", x2="1495", y2="721", stroke="#d2dcda")
    text(root, "ANCIENT ARCHAEOLOGICAL PANEL", panel_x, 763, 16, ACCENT, "bold")
    text(root, str(ancient_count), panel_x, 824, 50, ACCENT, "bold")
    text(root, "individuals · four site markers", panel_x, 856, 20)
    for index, line in enumerate(
        [
            "Liangdao: Larena et al. (2021)",
            "Suogang: Yang et al. (2020)",
            "Hanben / Gongguan: Wang et al. (2021)",
        ]
    ):
        text(root, line, panel_x, 901 + index * 29, 16)
    for index, line in enumerate(
        [
            "Counts describe individuals reported",
            "in the archaeological study panel,",
            "not the number of complete public",
            "mitochondrial consensus sequences.",
            "",
            "Markers show archaeological sites,",
            "not individual burial GPS positions.",
            "Suogang is rounded to 0.1 degree.",
        ]
    ):
        text(root, line, panel_x, 1015 + index * 25, 16, MUTED)
    text(root, "Made with Natural Earth · 1:50m land · public domain", 80, 1322, 17, MUTED)
    text(
        root,
        "Site sources: Taiwan Ministry of Culture heritage record; Yang 2020 Table 1; Wang 2021 sample metadata.",
        80,
        1352,
        16,
        MUTED,
    )
    text(
        root,
        "Liangdao's published rounded study coordinates (26.3 N, 120.2 E) are retained alongside the official site location in data/sites.csv.",
        80,
        1381,
        15,
        MUTED,
    )
    OUTPUT.mkdir(parents=True, exist_ok=True)
    svg_path = OUTPUT / "sample-sites.svg"
    png_path = OUTPUT / "sample-sites.png"
    ET.indent(root, space="  ")
    ET.ElementTree(root).write(svg_path, encoding="utf-8", xml_declaration=True)
    import cairosvg

    cairosvg.svg2png(url=str(svg_path), write_to=str(png_path), output_width=1550)
    print(
        json.dumps(
            {
                "svg": str(svg_path),
                "png": str(png_path),
                "ancient_individuals": ancient_count,
                "modern_individuals": 550,
            }
        )
    )


if __name__ == "__main__":
    main()
