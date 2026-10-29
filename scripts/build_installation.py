"""Build the exact-data Stump of Life installation as a metre-scale Y-up GLB.

Run: uv run --group geometry --group preview python scripts/build_installation.py

Scientific inputs stay separate from artistic geometry: the cap holes come from
the current holes.csv; six pending ancient tracks are dashed and unperforated.
The Liangdao-centred map is a thematic installation anchor, not an approved
exhibition location. This prototype establishes no structural/load rating.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import trimesh
import typer
from PIL import Image, ImageDraw, ImageFont
from shapely import geometry as sg
from shapely.ops import unary_union

from stump_of_life.io import read_fasta

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts"
PLINTH_HEIGHT = 0.025
WOOD_HEIGHT = 0.450
CAP_GAP = 0.008
CAP_THICKNESS = 0.003
CAP_BOTTOM = PLINTH_HEIGHT + WOOD_HEIGHT + CAP_GAP
CAP_TOP = CAP_BOTTOM + CAP_THICKNESS
WEST, EAST, SOUTH, NORTH = 117.5, 123.5, 21.0, 27.0
LON_SCALE = 3.0 / (EAST - WEST)
LAT_SCALE = LON_SCALE / math.cos(math.radians(24.0))
PENDING = [
    (55.0, "Liangdao 2"),
    (92.0, "Suogang B1"),
    (100.0, "Suogang B3"),
    (125.0, "Gongguan I13721"),
    (145.0, "Hanben I3620"),
    (160.0, "Hanben I8081"),
]
Y_UP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]], dtype=float)


def pbr(name: str, color: list[float], metal: float = 0.0, roughness: float = 0.85) -> Any:
    return trimesh.visual.material.PBRMaterial(
        name=name,
        baseColorFactor=color,
        metallicFactor=metal,
        roughnessFactor=roughness,
        doubleSided=False,
    )


MATERIALS = {
    "Map_Platform": pbr("Muted blue grey map plinth", [0.52, 0.64, 0.65, 1]),
    "Map_Land": pbr("Warm pale geographical land relief", [0.63, 0.65, 0.53, 1]),
    "Wood_Body": pbr("Natural bark with vertex colour", [1.0, 1.0, 1.0, 1.0], roughness=0.96),
    "Steel_Cap": pbr("Brushed silver prototype metal", [0.61, 0.66, 0.68, 1], 0.95, 0.29),
    "Track_Engravings": pbr("Dark engraved track marks", [0.10, 0.15, 0.15, 1], 0.55, 0.65),
    "Copper_Roots": pbr("Dark copper provenance roots", [0.40, 0.22, 0.11, 1], 0.55, 0.62),
    "Cap_Supports": pbr("Dark removable cap stand-offs", [0.14, 0.16, 0.15, 1], 0.65, 0.55),
}


def project(longitude: float, latitude: float, anchor: dict[str, str]) -> tuple[float, float]:
    """Source construction axes are X east, Y north, Z up, in metres."""
    return (
        (longitude - float(anchor["longitude"])) * LON_SCALE,
        (latitude - float(anchor["latitude"])) * LAT_SCALE,
    )


def polygons(shape: Any) -> list[Any]:
    if shape.is_empty:
        return []
    if shape.geom_type == "Polygon":
        return [shape]
    return [item for part in shape.geoms for item in polygons(part)]


def extrude(shape: Any, height: float, z: float) -> list[trimesh.Trimesh]:
    result = []
    for polygon in polygons(shape):
        mesh = trimesh.creation.extrude_polygon(polygon, height=height, engine="earcut")
        mesh.apply_translation([0, 0, z])
        result.append(mesh)
    return result


def flat(shape: Any, z: float) -> list[trimesh.Trimesh]:
    result = []
    for polygon in polygons(shape):
        vertices, faces = trimesh.creation.triangulate_polygon(polygon, engine="earcut")
        xyz = np.column_stack([vertices, np.full(len(vertices), z)])
        mesh = trimesh.Trimesh(vertices=xyz, faces=faces, process=False)
        if mesh.face_normals[0, 2] < 0:
            mesh.invert()
        result.append(mesh)
    return result


def annulus(radius: float, width: float = 0.00024) -> Any:
    return (
        sg.Point(0, 0)
        .buffer(radius + width / 2, quad_segs=128)
        .difference(sg.Point(0, 0).buffer(radius - width / 2, quad_segs=128))
    )


def dashed_track(radius: float) -> list[Any]:
    shapes = []
    for index in range(64):
        start = 2 * math.pi * index / 64
        end = start + 2 * math.pi * 0.56 / 64
        angles = np.linspace(start, end, 5)
        outer = [
            (math.cos(a) * (radius + 0.00018), math.sin(a) * (radius + 0.00018)) for a in angles
        ]
        inner = [
            (math.cos(a) * (radius - 0.00018), math.sin(a) * (radius - 0.00018))
            for a in angles[::-1]
        ]
        shapes.append(sg.Polygon(outer + inner))
    return shapes


def bark_body() -> trimesh.Trimesh:
    segments, levels = 128, 13
    vertices: list[list[float]] = []
    colors: list[list[int]] = []
    for level in range(levels):
        t = level / (levels - 1)
        for index in range(segments):
            angle = 2 * math.pi * index / segments
            radius = (
                0.269
                + 0.005 * math.sin(7 * angle)
                + 0.003 * math.sin(19 * angle + 0.4)
                + 0.0025 * math.sin(43 * angle)
                + 0.002 * math.sin(3 * angle + t)
            )
            z = PLINTH_HEIGHT + WOOD_HEIGHT * t + 0.0008 * t * math.sin(11 * angle)
            vertices.append([radius * math.cos(angle), radius * math.sin(angle), z])
            grain = 0.76 + 0.20 * math.sin(17 * angle + 0.2 * math.sin(t * 9))
            grain += 0.08 * math.sin(39 * angle + t * 3) + 0.04 * math.sin(t * 37 + angle)
            colors.append([int(145 * grain), int(101 * grain), int(60 * grain), 255])
    faces: list[list[int]] = []
    for level in range(levels - 1):
        for index in range(segments):
            following = (index + 1) % segments
            a, b = level * segments + index, level * segments + following
            c, d = a + segments, b + segments
            faces.extend([[a, b, d], [a, d, c]])
    bottom, top = len(vertices), len(vertices) + 1
    vertices.extend([[0, 0, PLINTH_HEIGHT], [0, 0, PLINTH_HEIGHT + WOOD_HEIGHT]])
    colors.extend([[116, 82, 48, 255], [139, 108, 65, 255]])
    for index in range(segments):
        following = (index + 1) % segments
        faces.append([bottom, following, index])
        a, b = (levels - 1) * segments + index, (levels - 1) * segments + following
        faces.append([top, a, b])
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    visual = trimesh.visual.TextureVisuals(material=MATERIALS["Wood_Body"])
    visual.vertex_attributes["color"] = np.asarray(colors, dtype=np.uint8)
    mesh.visual = visual
    return mesh


def root_tube(
    end: tuple[float, float],
    count: int,
    name: str,
    control_override: tuple[tuple[float, float], tuple[float, float]] | None = None,
) -> trimesh.Trimesh:
    control = {
        "Hanben": ((0.35, -0.30), (0.60, -0.65)),
        "Gongguan": ((0.50, -0.10), (1.25, -1.45)),
        "Suogang": ((-0.55, -0.35), (-0.70, -1.05)),
    }
    first, second = control[name] if control_override is None else control_override
    control_points = np.array([[0, 0], first, second, end], dtype=float)
    positions = []
    for t in np.linspace(0, 1, 73):
        positions.append(
            (1 - t) ** 3 * control_points[0]
            + 3 * (1 - t) ** 2 * t * control_points[1]
            + 3 * (1 - t) * t * t * control_points[2]
            + t**3 * control_points[3]
        )
    positions = np.asarray(positions)
    tangents = np.gradient(positions, axis=0)
    vertices = []
    half_width = 0.003 + 0.0008 * math.sqrt(count)
    for index, (point, tangent) in enumerate(zip(positions, tangents, strict=True)):
        perpendicular = np.array([-tangent[1], tangent[0]])
        perpendicular /= np.linalg.norm(perpendicular)
        taper = 1.0 - 0.58 * (index / (len(positions) - 1))
        for section in range(8):
            angle = 2 * math.pi * section / 8
            xy = point + perpendicular * half_width * taper * math.cos(angle)
            vertices.append([xy[0], xy[1], 0.0275 + 0.0015 * math.sin(angle)])
    faces = []
    for index in range(len(positions) - 1):
        for section in range(8):
            following = (section + 1) % 8
            a, b, c, d = (
                index * 8 + section,
                index * 8 + following,
                (index + 1) * 8 + section,
                (index + 1) * 8 + following,
            )
            faces.extend([[a, b, d], [a, d, c]])
    faces.extend([[0, index + 1, index] for index in range(1, 7)])
    last = (len(positions) - 1) * 8
    faces.extend([[last, last + index, last + index + 1] for index in range(1, 7)])
    return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)


def label_atlas(labels: list[str]) -> tuple[Image.Image, list[tuple[float, float]]]:
    atlas = Image.new("RGBA", (1536, 128 * len(labels)), (0, 0, 0, 0))
    draw = ImageDraw.Draw(atlas)
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 50)
    for index, label in enumerate(labels):
        size = min(50, int(50 * 1496 / max(font.getlength(label), 1)))
        fitted = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", size)
        draw.text((20, index * 128 + 25), label, font=fitted, fill=(35, 53, 55, 255))
    return atlas, [(index / len(labels), (index + 1) / len(labels)) for index in range(len(labels))]


def label_quad(
    center: tuple[float, float],
    width: float,
    height: float,
    uv_band: tuple[float, float],
    material: Any,
) -> trimesh.Trimesh:
    x, y = center
    vertices = [
        [x - width / 2, y + height / 2, 0.0297],
        [x + width / 2, y + height / 2, 0.0297],
        [x + width / 2, y - height / 2, 0.0297],
        [x - width / 2, y - height / 2, 0.0297],
    ]
    mesh = trimesh.Trimesh(vertices=vertices, faces=[[0, 2, 1], [0, 3, 2]], process=False)
    lo, hi = uv_band
    mesh.visual = trimesh.visual.TextureVisuals(
        uv=np.array([[0, 1 - lo], [1, 1 - lo], [1, 1 - hi], [0, 1 - hi]]), material=material
    )
    return mesh


def set_vertex_normals(mesh: trimesh.Trimesh) -> None:
    """Compute angle-weighted normals with NumPy, without optional SciPy."""
    triangle = mesh.vertices[mesh.faces]
    face_normals = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    lengths = np.linalg.norm(face_normals, axis=1)
    face_normals /= np.maximum(lengths[:, None], 1e-15)
    normals = np.zeros_like(mesh.vertices)
    for index in range(3):
        first = triangle[:, (index + 1) % 3] - triangle[:, index]
        second = triangle[:, (index + 2) % 3] - triangle[:, index]
        denominator = np.linalg.norm(first, axis=1) * np.linalg.norm(second, axis=1)
        cosine = np.sum(first * second, axis=1) / np.maximum(denominator, 1e-15)
        angle = np.arccos(np.clip(cosine, -1, 1))
        np.add.at(normals, mesh.faces[:, index], face_normals * angle[:, None])
    normals /= np.maximum(np.linalg.norm(normals, axis=1)[:, None], 1e-15)
    mesh.vertex_normals = normals


def hard_horizontal_edges(mesh: trimesh.Trimesh) -> trimesh.Trimesh:
    """Preserve positions while separating planar tops from edge shading."""
    normals = mesh.face_normals[:, 2]
    surfaces = []
    for mask in (normals > 0.99, normals < -0.99, np.abs(normals) <= 0.99):
        if not np.any(mask):
            continue
        surface = trimesh.Trimesh(
            vertices=mesh.vertices.copy(), faces=mesh.faces[mask], process=False
        )
        surface.remove_unreferenced_vertices()
        surfaces.append(surface)
    return trimesh.util.concatenate(surfaces)


def digital_dna_texture(geometry: dict[str, Any]) -> tuple[Image.Image, dict[str, Any]]:
    """Typeset actual MSA calls on a readable digital cap; never synthesize ancient calls."""
    aligned_path = ROOT / "output" / "taiwan" / "aligned.fasta"
    records = read_fasta(aligned_path, aligned=True)
    reference_id = geometry["reference_id"]
    # The metal export reports human-readable one-based MSA columns.
    columns = [column - 1 for column in geometry["angular_alignment_columns"]]
    if len(columns) != 175:
        raise ValueError("The reviewed Livistone texture expects 175 shared variable columns.")
    image = Image.new("RGB", (4096, 4096), (211, 220, 219))
    draw = ImageDraw.Draw(image)
    centre, scale = 2048, 7.6
    font_path = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
    small = ImageFont.truetype(font_path, 22)
    bold = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 34)
    reference_font = ImageFont.truetype(font_path, 11)
    colours = {"A": (0, 91, 66), "C": (19, 65, 160), "G": (131, 79, 0), "T": (162, 24, 67)}
    counts = {"substitution": 0, "insertion": 0, "deletion": 0, "unknown": 0}
    tracks = [(30.0, reference_id)] + [
        (float(row["radius_mm"]), row["sample_id"]) for row in geometry["guides"][1:]
    ]
    for radius_mm, sample_id in tracks:
        radius = radius_mm * scale
        draw.ellipse(
            (centre - radius, centre - radius, centre + radius, centre + radius),
            outline=(144, 163, 160),
            width=2,
        )
        for index, column in enumerate(columns):
            angle = 2 * math.pi * index / len(columns)
            x, y = centre + radius * math.sin(angle), centre - radius * math.cos(angle)
            ref, value = records[reference_id][column], records[sample_id][column]
            colour, font, glyph = (95, 111, 110), small, value
            if sample_id == reference_id:
                font = reference_font
            elif ref in "ACGT" and value in "ACGT" and ref != value:
                counts["substitution"] += 1
                colour, font = colours[value], bold
                draw.ellipse(
                    (x - 23, y - 23, x + 23, y + 23), fill=(235, 239, 235), outline=colour, width=2
                )
            elif ref == "-" and value in "ACGT":
                counts["insertion"] += 1
                colour = colours[value]
                draw.polygon(
                    [(x, y - 24), (x + 24, y), (x, y + 24), (x - 24, y)], outline=colour, width=2
                )
            elif value == "-" and ref in "ACGT":
                counts["deletion"] += 1
                glyph, font = "×", bold
            elif ref not in "ACGT-" or value not in "ACGT-":
                counts["unknown"] += 1
                colour = (135, 141, 140)
            draw.text((x, y), glyph, font=font, fill=colour, anchor="mm")
    for radius_mm, name in PENDING:
        radius = radius_mm * scale
        for dash in range(72):
            draw.arc(
                (centre - radius, centre - radius, centre + radius, centre + radius),
                start=dash * 5,
                end=dash * 5 + 2.8,
                fill=(105, 130, 122),
                width=3,
            )
        label_font = ImageFont.truetype(font_path, 25)
        status = "mtDNA unavailable" if name.startswith("Suogang ") else "sequence pending"
        label = f"{name}: {status}"
        box = draw.textbbox((0, 0), label, font=label_font)
        label_y = centre + radius * 0.68
        draw.rectangle(
            (centre - (box[2] + 20) / 2, label_y - 19, centre + (box[2] + 20) / 2, label_y + 20),
            fill=(211, 220, 219),
        )
        draw.text((centre, label_y), label, font=label_font, fill=(64, 91, 81), anchor="mm")
    for text, y, size in [
        ("LIANGDAO 1", centre - 60, 57),
        ("KF540505.1", centre + 5, 39),
        ("175 shared variable columns", centre + 60, 25),
        ("8,060–8,320 cal BP", centre + 100, 25),
    ]:
        font = ImageFont.truetype(font_path, size)
        draw.text((centre, y), text, font=font, fill=(35, 69, 60), anchor="mm")
    if counts["substitution"] != 294:
        raise ValueError(f"Texture call count disagrees with the reviewed 294 SNPs: {counts}.")
    return image, {
        "reference_id": reference_id,
        "shared_alignment_columns_zero_based": columns,
        "source_alignment_columns_one_based": geometry["angular_alignment_columns"],
        "coordinate_mode": "175 compressed shared variable MSA columns; includes insertion columns",
        "profile_tracks": [
            {"sample_id": key, "radius_mm_in_physical_plan": radius} for radius, key in tracks
        ],
        "pending_tracks": [
            {
                "sample": name,
                "radius_mm_in_physical_plan": radius,
                "status": "mtDNA unavailable; no MT reads in public deposit"
                if name.startswith("Suogang ")
                else "sequence pending",
            }
            for radius, name in PENDING
        ],
        "calls_relative_to_reference": counts,
        "aligned_fasta_sha256": hashlib.sha256(aligned_path.read_bytes()).hexdigest(),
        "texture_dimensions": [4096, 4096],
        "ancient_missing_calls_synthesized": False,
    }


def build_livistone() -> None:
    """Export an enlarged river island edition, preserving rigid geographic positions."""
    with (ROOT / "data" / "sites.csv").open(newline="") as handle:
        sites = list(csv.DictReader(handle))
    raw_geography = (ROOT / "data" / "geography" / "ne_50m_land.geojson").read_bytes()
    earth = json.loads(raw_geography)
    land_polygons = [
        polygon
        for feature in earth["features"]
        for polygon in polygons(sg.shape(feature["geometry"]))
    ]
    taiwan = next(polygon for polygon in land_polygons if polygon.contains(sg.Point(121, 23.7)))
    cosine = math.cos(math.radians(24))
    points = np.asarray(taiwan.exterior.coords).copy()
    points[:, 0] *= cosine
    origin = points.mean(axis=0)
    angle = math.radians(83.75)
    rotation = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
    rotated = (points - origin) @ rotation
    metres_per_degree = 18.5 / np.ptp(rotated[:, 0])
    rotated *= metres_per_degree
    sites_raw = np.array(
        [[float(row["longitude"]) * cosine, float(row["latitude"])] for row in sites]
    )
    site_positions = (sites_raw - origin) @ rotation * metres_per_degree
    radii = {"Liangdao": 2.25, "Suogang": 0.8, "Hanben": 0.0, "Gongguan": 0.7}
    low, high = rotated.min(axis=0), rotated.max(axis=0)
    for row, position in zip(sites, site_positions, strict=True):
        radius = radii[row["site"]]
        low, high = np.minimum(low, position - radius), np.maximum(high, position + radius)
    centre = (low + high) / 2
    rotated -= centre
    site_positions -= centre
    # Construction XY is a rigidly rotated East/North plane; Y-up GLB uses Z=-Y.
    main_polygon = sg.Polygon(rotated)
    positions = {row["site"]: position for row, position in zip(sites, site_positions, strict=True)}
    matsu = positions["Liangdao"]
    shapes = {"Main_Taiwan_Island": main_polygon}
    name_to_site = {
        "Matsu_Islet": "Liangdao",
        "Penghu_Islet": "Suogang",
        "Green_Island_Islet": "Gongguan",
    }
    for name, site in name_to_site.items():
        point, radius = positions[site], radii[site]
        boundary = []
        for index in range(64):
            a = 2 * math.pi * index / 64
            r = radius * (0.96 + 0.025 * math.sin(5 * a) + 0.012 * math.sin(9 * a + 0.4))
            boundary.append([point[0] + r * math.cos(a), point[1] + r * math.sin(a)])
        shapes[name] = sg.Polygon(boundary)
    earth_material = pbr("River island sand and stone", [0.43, 0.37, 0.26, 1.0], 0.0, 0.98)
    grass_material = pbr("Walkable Taiwan island greenery", [0.27, 0.43, 0.23, 1.0], 0.0, 0.96)
    copper_material = pbr(
        "Copper geographical provenance roots", [0.38, 0.18, 0.07, 1.0], 0.55, 0.55
    )
    scene = trimesh.Scene()
    boundaries = {}
    for name, shape in shapes.items():
        mesh = trimesh.util.concatenate(extrude(shape, 2.75, -2.5))
        mesh = hard_horizontal_edges(mesh)
        mesh.visual = trimesh.visual.TextureVisuals(material=earth_material)
        mesh.apply_transform(Y_UP)
        set_vertex_normals(mesh)
        scene.add_geometry(mesh, geom_name=name, node_name=name)
        boundaries[name] = [[float(x), float(-y)] for x, y in shape.exterior.coords]
        interior = shape.buffer(-0.20 if name == "Main_Taiwan_Island" else -0.14)
        grass = trimesh.util.concatenate(flat(interior, 0.251))
        grass.visual = trimesh.visual.TextureVisuals(material=grass_material)
        grass.apply_transform(Y_UP)
        set_vertex_normals(grass)
        scene.add_geometry(grass, geom_name=f"{name}_Grass", node_name=f"{name}_Grass")
    wood = bark_body()
    transform = np.eye(4)
    transform[0, 0] = transform[1, 1] = 1.75 / 0.269
    transform[2, 2] = 0.630 / WOOD_HEIGHT
    transform[:3, 3] = [matsu[0], matsu[1], 0.25 - PLINTH_HEIGHT * transform[2, 2]]
    wood.apply_transform(transform)
    wood.apply_transform(Y_UP)
    set_vertex_normals(wood)
    scene.add_geometry(wood, geom_name="Giant_Wood_Stump", node_name="Giant_Wood_Stump")
    cap_radius, cap_bottom, cap_top = 1.70, 0.915, 0.94
    cap_shape = sg.Point(float(matsu[0]), float(matsu[1])).buffer(cap_radius, quad_segs=64)
    closed_cap = trimesh.util.concatenate(extrude(cap_shape, cap_top - cap_bottom, cap_bottom))
    source_watertight = closed_cap.is_watertight
    normal_z = closed_cap.face_normals[:, 2]
    top_mesh = trimesh.Trimesh(
        vertices=closed_cap.vertices.copy(), faces=closed_cap.faces[normal_z > 0.99], process=False
    )
    top_mesh.remove_unreferenced_vertices()
    cap_body = trimesh.Trimesh(
        vertices=closed_cap.vertices.copy(), faces=closed_cap.faces[normal_z <= 0.99], process=False
    )
    cap_body.remove_unreferenced_vertices()
    cap_body = hard_horizontal_edges(cap_body)
    cap_body.visual = trimesh.visual.TextureVisuals(material=MATERIALS["Steel_Cap"])
    geometry = json.loads((ROOT / "output" / "taiwan" / "metal" / "geometry.json").read_text())
    texture, texture_metadata = digital_dna_texture(geometry)
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    texture_path = ARTIFACTS / "stump-of-life-livistone-cap.png"
    texture.save(texture_path, optimize=True)
    texture_material = trimesh.visual.material.PBRMaterial(
        name="Real mitochondrial DNA calls; solid engraved digital cap",
        baseColorFactor=[1.0, 1.0, 1.0, 1.0],
        baseColorTexture=texture,
        metallicFactor=0.10,
        roughnessFactor=0.82,
        doubleSided=False,
    )
    top_mesh.visual = trimesh.visual.TextureVisuals(
        uv=np.column_stack(
            [
                (top_mesh.vertices[:, 0] - matsu[0]) / (2 * cap_radius) + 0.5,
                (top_mesh.vertices[:, 1] - matsu[1]) / (2 * cap_radius) + 0.5,
            ]
        ),
        material=texture_material,
    )
    for name, mesh in [("Giant_Solid_Cap_Edge", cap_body), ("Giant_Solid_DNA_Cap", top_mesh)]:
        mesh.apply_transform(Y_UP)
        set_vertex_normals(mesh)
        scene.add_geometry(mesh, geom_name=name, node_name=name)
    # Curves link the exact rigidly transformed endpoints; their controls remain
    # inside the planned land/islet envelope rather than stretching the channel fit.
    scale_ratio = metres_per_degree / (LON_SCALE / cosine)
    root_meshes = []
    for row in sites:
        if row["site"] == "Liangdao":
            continue
        delta = positions[row["site"]] - matsu
        perpendicular = np.array([-delta[1], delta[0]]) / np.linalg.norm(delta)
        controls = [
            matsu + delta * 0.3 + perpendicular * 0.55,
            matsu + delta * 0.7 + perpendicular * 0.28,
        ]
        controls = [np.clip(point, low - centre + 0.14, high - centre - 0.14) for point in controls]
        relative_controls = tuple(tuple((point - matsu) / scale_ratio) for point in controls)
        root = root_tube(
            tuple(delta / scale_ratio),
            int(row["ancient_individuals"]),
            row["site"],
            control_override=relative_controls,
        )
        mapping = np.eye(4)
        mapping[:2, :2] *= scale_ratio
        mapping[0, 3], mapping[1, 3] = matsu
        mapping[2, 3] = 0.2475
        root.apply_transform(mapping)
        root_meshes.append(root)
    roots = trimesh.util.concatenate(root_meshes)
    roots.visual = trimesh.visual.TextureVisuals(material=copper_material)
    roots.apply_transform(Y_UP)
    set_vertex_normals(roots)
    scene.add_geometry(
        roots,
        geom_name="Copper_Provenance_Roots_Not_Bridges",
        node_name="Copper_Provenance_Roots_Not_Bridges",
    )
    labels = [
        "TAIWAN",
        "Liangdao / Matsu",
        "Suogang / Penghu",
        "Gongguan / Green Island",
        "Hanben / Blihun: 45 individuals",
    ]
    atlas, bands = label_atlas(labels)
    labels_material = trimesh.visual.material.PBRMaterial(
        name="Digital island site labels",
        baseColorFactor=[1.0, 1.0, 1.0, 1.0],
        baseColorTexture=atlas,
        metallicFactor=0,
        roughnessFactor=1,
        alphaMode="MASK",
        alphaCutoff=0.3,
        doubleSided=True,
    )
    representative = main_polygon.representative_point()
    label_specs = [
        ((representative.x, representative.y), 5.0, 0.50),
        ((matsu[0] - 0.3, matsu[1] + 1.95), 3.0, 0.30),
        ((positions["Suogang"][0], positions["Suogang"][1]), 1.45, 0.21),
        ((positions["Gongguan"][0], positions["Gongguan"][1]), 1.30, 0.21),
        ((positions["Hanben"][0] - 0.4, positions["Hanben"][1] + 0.4), 3.6, 0.32),
    ]
    label_meshes = []
    for (point, width, height), band in zip(label_specs, bands, strict=True):
        label = label_quad(point, width, height, band, labels_material)
        label.apply_translation([0, 0, 0.254 - 0.0297])
        label_meshes.append(label)
    labels_mesh = trimesh.util.concatenate(label_meshes)
    labels_mesh.apply_transform(Y_UP)
    set_vertex_normals(labels_mesh)
    scene.add_geometry(
        labels_mesh, geom_name="Geographical_Island_Labels", node_name="Geographical_Island_Labels"
    )
    from shapely.ops import nearest_points

    boundary = nearest_points(main_polygon.boundary, sg.Point(matsu))[0]
    boundary_point = np.array(boundary.coords[0])
    inside = np.array(main_polygon.representative_point().coords[0])
    inside = (
        boundary_point + (inside - boundary_point) / np.linalg.norm(inside - boundary_point) * 0.5
    )
    approach = matsu + (boundary_point - matsu) / np.linalg.norm(boundary_point - matsu) * (
        cap_radius - 0.34
    )
    final_path = ARTIFACTS / "stump-of-life-livistone.glb"
    scene.export(str(final_path), file_type="glb", include_normals=True)
    metadata = {
        "asset": final_path.name,
        "edition": "Livistone enlarged river island; separate from the physical seat",
        "units": "metres",
        "up_axis": "Y",
        "island_bottom_y": -2.5,
        "island_top_y": 0.25,
        "bounds_min": scene.bounds[0].tolist(),
        "bounds_max": scene.bounds[1].tolist(),
        "footprint_is_centred": True,
        "planned_footprint_xz_m": [float(high[0] - low[0]), float(high[1] - low[1])],
        "geographic_transform": {
            "projection": "Equirectangular at24N",
            "longitude_cosine": cosine,
            "projected_origin": origin.tolist(),
            "rotation_degrees": 83.75,
            "uniform_metres_per_projected_degree": metres_per_degree,
            "projected_centre_offset": centre.tolist(),
            "y_up_convention": "localX=rotatedEastNorthX; localZ=-rotatedEastNorthY",
            "geography_warped": False,
            "main_taiwan_length_m": 18.5,
        },
        "land_boundaries_xz": boundaries,
        "walkable_meshes": list(shapes),
        "site_anchors": [
            {
                "site": row["site"],
                "latitude": float(row["latitude"]),
                "longitude": float(row["longitude"]),
                "individuals": int(row["ancient_individuals"]),
                "local_xz": [float(point[0]), float(-point[1])],
                "symbolically_enlarged_islet_radius_m": radii[row["site"]],
            }
            for row, point in zip(sites, site_positions, strict=True)
        ],
        "stump": {
            "centre_xz": [float(matsu[0]), float(-matsu[1])],
            "body_radius_m": 1.75,
            "cap_radius_m": cap_radius,
            "cap_top_y": cap_top,
            "cap_bottom_y": cap_bottom,
            "solid_cap_source_watertight": source_watertight,
            "through_holes": 0,
            "meaning": "Liangdao/Matsu reference anchor, not an exhibition venue GPS claim",
        },
        "bridge_targets": {
            "main_coast_boundary": [float(boundary_point[0]), 0.25, float(-boundary_point[1])],
            "main_walkable_approach": [float(inside[0]), 0.25, float(-inside[1])],
            "cap_walkable_approach": [float(approach[0]), cap_top, float(-approach[1])],
            "copper_roots_are_not_walkway_bridges": True,
        },
        "dna_texture": {
            **texture_metadata,
            "file": texture_path.name,
            "sha256": hashlib.sha256(texture_path.read_bytes()).hexdigest(),
        },
        "natural_earth_sha256": hashlib.sha256(raw_geography).hexdigest(),
        "modern_collection_points_invented": False,
        "vertex_count": sum(len(m.vertices) for m in scene.geometry.values()),
        "triangle_count": sum(len(m.faces) for m in scene.geometry.values()),
        "file_bytes": final_path.stat().st_size,
        "interpretation": "Root curves are geographic provenance artwork, not migration routes or phylogeny",
    }
    (ARTIFACTS / "stump-of-life-livistone.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n"
    )
    print(
        json.dumps(
            {
                key: metadata[key]
                for key in [
                    "asset",
                    "bounds_min",
                    "bounds_max",
                    "vertex_count",
                    "file_bytes",
                    "stump",
                    "bridge_targets",
                ]
            },
            indent=2,
        )
    )


def main(
    cap_treatment: str = "holes", engraved_mark_diameter_mm: float = 1.5, edition: str = "physical"
) -> None:
    if edition == "livistone":
        build_livistone()
        return
    if edition != "physical":
        raise ValueError("edition must be 'physical' or 'livistone'.")
    if cap_treatment not in {"holes", "engraved"}:
        raise ValueError("cap_treatment must be 'holes' or 'engraved'.")
    if engraved_mark_diameter_mm <= 0:
        raise ValueError("engraved_mark_diameter_mm must be positive.")
    with (ROOT / "data" / "sites.csv").open(newline="") as handle:
        sites = list(csv.DictReader(handle))
    anchor = next(site for site in sites if site["site"] == "Liangdao")
    holes_path = ROOT / "output" / "taiwan" / "metal" / "holes.csv"
    with holes_path.open(newline="") as handle:
        holes = list(csv.DictReader(handle))
    geometry = json.loads((ROOT / "output" / "taiwan" / "metal" / "geometry.json").read_text())
    if len(holes) != 294 or geometry["hole_count"] != 294:
        raise ValueError("This reviewed installation expects the current exact 294-hole design.")
    parts: dict[str, list[trimesh.Trimesh]] = defaultdict(list)
    hole_shapes = [
        sg.Point((float(row["x_mm"]) - 250) / 1000, -(float(row["y_mm"]) - 250) / 1000).buffer(
            float(row["hole_diameter_mm"]) / 2000, quad_segs=4
        )
        for row in holes
    ]
    perforated = sg.Point(0, 0).buffer(0.25, quad_segs=128).difference(unary_union(hole_shapes))
    if perforated.geom_type != "Polygon" or len(perforated.interiors) != 294:
        raise ValueError("Cap Boolean did not preserve 294 distinct through-holes.")
    cap_shape = (
        perforated if cap_treatment == "holes" else sg.Point(0, 0).buffer(0.25, quad_segs=128)
    )
    parts["Steel_Cap"].extend(extrude(cap_shape, CAP_THICKNESS, CAP_BOTTOM))
    track_shapes = [annulus(0.030)]
    modern_radii = sorted({float(row["radius_mm"]) for row in holes})
    track_shapes.extend(annulus(radius / 1000) for radius in modern_radii)
    for radius, _ in PENDING:
        track_shapes.extend(dashed_track(radius / 1000))
    for track in track_shapes:
        parts["Track_Engravings"].extend(flat(track.intersection(cap_shape), CAP_TOP + 0.000018))
    if cap_treatment == "engraved":
        for row in holes:
            mark = sg.Point(
                (float(row["x_mm"]) - 250) / 1000,
                -(float(row["y_mm"]) - 250) / 1000,
            ).buffer(engraved_mark_diameter_mm / 2000, quad_segs=4)
            parts["Track_Engravings"].extend(flat(mark, CAP_TOP + 0.000019))
    for angle in np.linspace(0, 2 * math.pi, 4, endpoint=False):
        support = trimesh.creation.cylinder(radius=0.006, height=CAP_GAP, sections=12)
        support.apply_translation(
            [0.15 * math.cos(angle), 0.15 * math.sin(angle), CAP_BOTTOM - CAP_GAP / 2]
        )
        parts["Cap_Supports"].append(support)
    minimum_x, minimum_y = project(WEST, SOUTH, anchor)
    maximum_x, maximum_y = project(EAST, NORTH, anchor)
    frame = sg.box(minimum_x, minimum_y, maximum_x, maximum_y)
    rounded_frame = frame.buffer(-0.025).buffer(0.025, quad_segs=4)
    parts["Map_Platform"].extend(extrude(rounded_frame, PLINTH_HEIGHT, 0))
    raw_geography = (ROOT / "data" / "geography" / "ne_50m_land.geojson").read_bytes()
    earth = json.loads(raw_geography)
    geographic_clip = sg.box(WEST, SOUTH, EAST, NORTH)
    for feature in earth["features"]:
        clipped = sg.shape(feature["geometry"]).intersection(geographic_clip)
        for polygon in polygons(clipped):
            exterior = [project(lon, lat, anchor) for lon, lat, *_ in polygon.exterior.coords]
            interiors = [
                [project(lon, lat, anchor) for lon, lat, *_ in ring.coords]
                for ring in polygon.interiors
            ]
            parts["Map_Land"].extend(
                extrude(sg.Polygon(exterior, interiors), 0.0015, PLINTH_HEIGHT)
            )
    root_metadata = []
    for site in sites:
        x, y = project(float(site["longitude"]), float(site["latitude"]), anchor)
        count = int(site["ancient_individuals"])
        if site["site"] != "Liangdao":
            parts["Copper_Roots"].append(root_tube((x, y), count, site["site"]))
        marker = trimesh.creation.cylinder(radius=0.014, height=0.0016, sections=24)
        marker.apply_translation([x, y, 0.0275])
        parts["Copper_Roots"].append(marker)
        root_metadata.append(
            {
                "site": site["site"],
                "latitude": float(site["latitude"]),
                "longitude": float(site["longitude"]),
                "individual_count": count,
                "endpoint_y_up_metres": [x, 0.0275, -y],
                "meaning": "local thematic anchor"
                if site["site"] == "Liangdao"
                else "geographical provenance; not migration or phylogeny",
            }
        )
    scene = trimesh.Scene()
    for name, meshes in parts.items():
        mesh = trimesh.util.concatenate(meshes)
        if name in {"Steel_Cap", "Map_Platform", "Map_Land"}:
            mesh = hard_horizontal_edges(mesh)
        material = MATERIALS[name]
        if name == "Track_Engravings" and cap_treatment == "engraved":
            material = pbr("Dark laser-engraved data marks", [0.008, 0.008, 0.008, 1.0], 0.0, 0.94)
        mesh.visual = trimesh.visual.TextureVisuals(material=material)
        mesh.apply_transform(Y_UP)
        scene.add_geometry(mesh, geom_name=name, node_name=name)
    wood = bark_body()
    wood.apply_transform(Y_UP)
    scene.add_geometry(wood, geom_name="Wood_Body", node_name="Wood_Body")
    labels = [
        "STUMP OF LIFE",
        "Liangdao / Matsu · 2 ancient individuals",
        "Suogang / Penghu · 2",
        "Hanben / Blihun · 45",
        "Gongguan / Green Island · 1",
        "9 selected modern sequences · 6 ancient tracks pending",
    ]
    atlas, bands = label_atlas(labels)
    label_material = trimesh.visual.material.PBRMaterial(
        name="Embedded geographical label atlas",
        baseColorFactor=[1.0, 1.0, 1.0, 1.0],
        baseColorTexture=atlas,
        metallicFactor=0,
        roughnessFactor=1,
        alphaMode="MASK",
        alphaCutoff=0.3,
        doubleSided=True,
    )
    specs = [
        ((0.136, -2.68), 1.7, 0.14),
        ((-0.70, 0.19), 0.95, 0.095),
        ((-0.70, -1.67), 0.73, 0.075),
        ((1.15, -0.99), 0.75, 0.075),
        ((1.10, -2.12), 0.90, 0.09),
        ((0.136, -2.84), 2.4, 0.13),
    ]
    label_meshes = []
    for (center, width, height), band in zip(specs, bands, strict=True):
        label_meshes.append(label_quad(center, width, height, band, label_material))
    label_mesh = trimesh.util.concatenate(label_meshes)
    label_mesh.apply_transform(Y_UP)
    scene.add_geometry(label_mesh, geom_name="Map_Labels", node_name="Map_Labels")
    vertex_count = sum(len(mesh.vertices) for mesh in scene.geometry.values())
    if vertex_count >= 100000:
        raise ValueError(f"Asset exceeds the 100,000-vertex budget: {vertex_count}.")
    for mesh in scene.geometry.values():
        set_vertex_normals(mesh)
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    stem = "stump-of-life" if cap_treatment == "holes" else "stump-of-life-engraved"
    path = ARTIFACTS / f"{stem}.glb"
    scene.export(str(path), file_type="glb", include_normals=True)
    tracks = [
        {"sample": "Liangdao 1 / KF540505.1", "radius_mm": 30, "status": "published reference"}
    ]
    tracks.extend(
        {"sample": name, "radius_mm": radius, "status": "pending; dashed; no fabricated holes"}
        for radius, name in PENDING
    )
    tracks.extend(
        {
            "sample": guide["sample_id"],
            "radius_mm": guide["radius_mm"],
            "status": "real modern sequence",
        }
        for guide in geometry["guides"][1:]
    )
    metadata = {
        "asset": path.name,
        "units": "metres",
        "up_axis": "Y",
        "ground_y": 0,
        "stump_axis_local_xz": [0, 0],
        "bounds_min": scene.bounds[0].tolist(),
        "bounds_max": scene.bounds[1].tolist(),
        "vertex_count": vertex_count,
        "triangle_count": sum(len(mesh.faces) for mesh in scene.geometry.values()),
        "geometry_nodes": list(scene.geometry),
        "file_bytes": path.stat().st_size,
        "cap": {
            "treatment": cap_treatment,
            "diameter_m": 0.5,
            "prototype_thickness_m": CAP_THICKNESS,
            "bottom_y": CAP_BOTTOM,
            "top_y": CAP_TOP,
            "distinct_through_holes": len(cap_shape.interiors),
            "engraved_substitution_marks": 294 if cap_treatment == "engraved" else 0,
            "engraved_mark_diameter_mm": engraved_mark_diameter_mm
            if cap_treatment == "engraved"
            else None,
            "visual_mark_depth": "Dark surface marks, not through-holes"
            if cap_treatment == "engraved"
            else "Actual through-holes",
            "source_solid_watertight": all(mesh.is_watertight for mesh in parts["Steel_Cap"]),
            "hole_circle_segments": 16,
            "source": str(holes_path.relative_to(ROOT)),
            "source_sha256": hashlib.sha256(holes_path.read_bytes()).hexdigest(),
            "coordinate_transform": "CSV x right/y down converted to construction x east/y north; then Z-up to Y-up",
            "tracks": tracks,
            "track_count": len(tracks),
            "pending_ancient_tracks": 6,
        },
        "body": {
            "nominal_diameter_m": 0.54,
            "wood_height_m": WOOD_HEIGHT,
            "cap_gap_m": CAP_GAP,
            "description": "illustrative natural bark geometry and vertex colours; removable cap concept",
        },
        "map": {
            "width_m": 3.0,
            "plinth_height_m": PLINTH_HEIGHT,
            "extent": [WEST, SOUTH, EAST, NORTH],
            "anchor_site": anchor,
            "thematic_anchor_not_exhibition_venue": True,
            "source": "Natural Earth 1:50m public-domain land",
            "source_sha256": hashlib.sha256(raw_geography).hexdigest(),
            "small_islands": "Site-coordinate markers are exact; islands smaller than source coastline resolution are not invented",
        },
        "root_endpoints": root_metadata,
        "modern_collection_coordinates_invented": False,
        "interpretation": "Geographical provenance artwork; root curves are designed shapes, not migration routes or phylogeny",
        "engineering_status": "Visual prototype; no load rating, structural grade or fabrication approval",
    }
    (ARTIFACTS / f"{stem}.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n"
    )
    print(
        json.dumps(
            {
                key: metadata[key]
                for key in (
                    "asset",
                    "bounds_min",
                    "bounds_max",
                    "vertex_count",
                    "triangle_count",
                    "file_bytes",
                )
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    typer.run(main)
