"""Render the exported installation GLB faithfully with Blender 5.

Run: blender --background --python scripts/render_installation.py
Portrait: blender --background --python scripts/render_installation.py -- --view cover
Only studio lighting, cameras, a ground plane and a procedural bark bump shader
are added. No scientific geometry or perforation is generated in this script.
The source GLB is a generated design prototype, not a selected real timber specimen.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
import sys
from pathlib import Path
from typing import Any

import bpy
from mathutils import Vector

ROOT = Path(__file__).resolve().parents[1]
ASSET = ROOT / "artifacts" / "stump-of-life.glb"
IMAGES = ROOT / "docs" / "images"
SCENES = ROOT / "output" / "installation"


def aim(obj: Any, target: tuple[float, float, float]) -> None:
    obj.rotation_euler = (Vector(target) - obj.location).to_track_quat("-Z", "Y").to_euler()


def light(
    name: str,
    position: tuple[float, float, float],
    target: tuple[float, float, float],
    energy: float,
    size: float,
    colour: tuple[float, float, float] = (1.0, 1.0, 1.0),
) -> None:
    data = bpy.data.lights.new(name, "AREA")
    data.energy = energy
    data.shape = "DISK"
    data.size = size
    data.color = colour
    obj = bpy.data.objects.new(name, data)
    bpy.context.collection.objects.link(obj)
    obj.location = position
    aim(obj, target)


def camera(
    name: str,
    position: tuple[float, float, float],
    target: tuple[float, float, float],
    scale: float,
) -> Any:
    data = bpy.data.cameras.new(name)
    data.type = "ORTHO"
    data.ortho_scale = scale
    data.clip_start = 0.01
    data.clip_end = 100
    obj = bpy.data.objects.new(name, data)
    bpy.context.collection.objects.link(obj)
    obj.location = position
    aim(obj, target)
    return obj


def bark_bump(darker: bool = False) -> None:
    for material in bpy.data.materials:
        if "Natural bark" not in material.name or not material.use_nodes:
            continue
        nodes = material.node_tree.nodes
        links = material.node_tree.links
        shader = nodes.get("Principled BSDF")
        if shader is None:
            continue
        coordinates = nodes.new("ShaderNodeTexCoord")
        stretch = nodes.new("ShaderNodeVectorMath")
        stretch.operation = "MULTIPLY"
        stretch.inputs[1].default_value = (3.0, 3.0, 0.22)
        noise = nodes.new("ShaderNodeTexNoise")
        noise.inputs["Scale"].default_value = 70.0
        noise.inputs["Detail"].default_value = 3.0
        bump = nodes.new("ShaderNodeBump")
        bump.inputs["Strength"].default_value = 0.36
        bump.inputs["Distance"].default_value = 0.0025
        links.new(coordinates.outputs["Generated"], stretch.inputs[0])
        links.new(stretch.outputs["Vector"], noise.inputs["Vector"])
        links.new(noise.outputs["Fac"], bump.inputs["Height"])
        links.new(bump.outputs["Normal"], shader.inputs["Normal"])
        # Material study only: preserve the asset's vertex-colour timber while
        # adding a restrained vertical bark grain, not a claimed timber photo.
        ramp = nodes.new("ShaderNodeValToRGB")
        ramp.color_ramp.elements[0].position = 0.28
        ramp.color_ramp.elements[0].color = (
            (0.18, 0.13, 0.08, 1) if darker else (0.28, 0.23, 0.18, 1)
        )
        ramp.color_ramp.elements[1].position = 0.74
        ramp.color_ramp.elements[1].color = (
            (0.68, 0.50, 0.30, 1) if darker else (1.0, 0.89, 0.70, 1)
        )
        multiply = nodes.new("ShaderNodeMixRGB")
        multiply.blend_type = "MULTIPLY"
        multiply.inputs[0].default_value = 0.70 if darker else 0.58
        original_links = list(shader.inputs["Base Color"].links)
        if original_links:
            links.new(original_links[0].from_socket, multiply.inputs[1])
        else:
            multiply.inputs[1].default_value = shader.inputs["Base Color"].default_value
        links.new(noise.outputs["Fac"], ramp.inputs["Fac"])
        links.new(ramp.outputs["Color"], multiply.inputs[2])
        links.new(multiply.outputs["Color"], shader.inputs["Base Color"])


def metal_grain(restrained: bool = False) -> None:
    """Render-only fine brushed grain, without altering data mark geometry."""
    for material in bpy.data.materials:
        if "Brushed silver" not in material.name or not material.use_nodes:
            continue
        nodes, links = material.node_tree.nodes, material.node_tree.links
        shader = nodes.get("Principled BSDF")
        if shader is None:
            continue
        coordinates = nodes.new("ShaderNodeTexCoord")
        stretch = nodes.new("ShaderNodeVectorMath")
        stretch.operation = "MULTIPLY"
        stretch.inputs[1].default_value = (1.0, 180.0, 1.0)
        noise = nodes.new("ShaderNodeTexNoise")
        noise.inputs["Scale"].default_value = 28.0
        noise.inputs["Detail"].default_value = 2.0
        bump = nodes.new("ShaderNodeBump")
        bump.inputs["Strength"].default_value = 0.05 if restrained else 0.10
        bump.inputs["Distance"].default_value = 0.000004 if restrained else 0.000014
        links.new(coordinates.outputs["Generated"], stretch.inputs[0])
        links.new(stretch.outputs["Vector"], noise.inputs["Vector"])
        links.new(noise.outputs["Fac"], bump.inputs["Height"])
        links.new(bump.outputs["Normal"], shader.inputs["Normal"])


def mesh_signature(meshes: list[Any]) -> str:
    """Fingerprint imported positions, face indices and UVs before/after staging."""
    digest = hashlib.sha256()
    for obj in sorted(meshes, key=lambda item: item.name):
        digest.update(obj.name.encode())
        for vertex in obj.data.vertices:
            digest.update(struct.pack("<3f", *vertex.co))
        for polygon in obj.data.polygons:
            digest.update(struct.pack(f"<{len(polygon.vertices)}I", *polygon.vertices))
        for layer in obj.data.uv_layers:
            digest.update(layer.name.encode())
            for point in layer.data:
                digest.update(struct.pack("<2f", *point.uv))
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cap-treatment", choices=("holes", "engraved"), default="holes")
    parser.add_argument("--edition", choices=("physical", "livistone"), default="physical")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--view", choices=("all", "cap", "top", "cover"), default="all")
    parser.add_argument(
        "--output",
        type=Path,
        help="PNG destination for --view cover; its Blender scene and provenance use the same stem.",
    )
    arguments = parser.parse_args(sys.argv[sys.argv.index("--") + 1 :] if "--" in sys.argv else [])
    digital = arguments.edition == "livistone"
    cover = arguments.view == "cover"
    if arguments.output is not None:
        if not cover:
            parser.error("--output applies to --view cover only.")
        if arguments.output.suffix.lower() != ".png":
            parser.error("--output must have a .png extension.")
    if cover:
        if digital:
            raise ValueError(
                "The portrait view uses the physical seat, not the enlarged digital edition."
            )
        arguments.cap_treatment = "engraved"
    suffix = (
        "-livistone" if digital else ("" if arguments.cap_treatment == "holes" else "-engraved")
    )
    asset = ROOT / "artifacts" / f"stump-of-life{suffix}.glb"
    if not asset.is_file():
        raise FileNotFoundError(asset)
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)
    bpy.context.preferences.filepaths.save_version = 0
    bpy.ops.import_scene.gltf(filepath=str(asset))
    meshes = [obj for obj in bpy.context.scene.objects if obj.type == "MESH"]
    original_signature = mesh_signature(meshes)
    original_asset_sha = hashlib.sha256(asset.read_bytes()).hexdigest()
    corners = [obj.matrix_world @ Vector(corner) for obj in meshes for corner in obj.bound_box]
    minimum = [min(corner[axis] for corner in corners) for axis in range(3)]
    maximum = [max(corner[axis] for corner in corners) for axis in range(3)]
    print(json.dumps({"imported_blender_z_up_bounds": [minimum, maximum]}), flush=True)
    if digital:
        expected_bounds = abs(minimum[2] + 2.5) < 0.001 and abs(maximum[2] - 0.94) < 0.001
    else:
        expected_bounds = 0.48 < maximum[2] < 0.50 and abs(minimum[2]) < 0.001
    if not expected_bounds:
        raise ValueError(
            "The imported GLB was not converted to the expected metre-scale Z-up scene."
        )
    scene = bpy.context.scene
    scene.unit_settings.system = "METRIC"
    scene.unit_settings.scale_length = 1.0
    scene.render.engine = "CYCLES"
    scene.cycles.device = "CPU"
    scene.cycles.samples = 8 if arguments.quick else (64 if cover else 40)
    scene.cycles.use_denoising = True
    scene.cycles.max_bounces = 6
    scene.render.threads_mode = "FIXED"
    scene.render.threads = 4
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGB"
    scene.render.image_settings.color_depth = "8"
    scene.render.film_transparent = False
    scene.view_settings.view_transform = "AgX"
    scene.view_settings.exposure = -0.15 if cover else -0.3
    world = bpy.data.worlds.new("Neutral exhibition studio")
    world.use_nodes = True
    world.node_tree.nodes["Background"].inputs["Color"].default_value = (
        (0.38, 0.43, 0.41, 1) if cover else (0.67, 0.73, 0.77, 1)
    )
    world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.12 if cover else 0.42
    scene.world = world
    bark_bump(darker=cover)
    metal_grain(restrained=cover)
    floor = bpy.data.materials.new("Neutral warm studio floor")
    floor.use_nodes = True
    shader = floor.node_tree.nodes.get("Principled BSDF")
    shader.inputs["Base Color"].default_value = (
        (0.08, 0.11, 0.10, 1)
        if cover
        else ((0.16, 0.27, 0.28, 1) if digital else (0.42, 0.40, 0.36, 1))
    )
    shader.inputs["Roughness"].default_value = 0.36 if digital else 0.94
    bpy.ops.mesh.primitive_plane_add(size=90, location=(0, -1.2, 0.0 if digital else -0.004))
    bpy.context.object.name = "Studio ground only; not part of exported installation"
    bpy.context.object.data.materials.append(floor)
    if digital:
        light("Digital island soft key", (-8, -10, 18), (0, 0, 0), 13000, 16)
        light("Digital island fill", (12, 8, 14), (2, 1, 0), 6500, 12)
        hero = camera("Enlarged Livistone river island", (23, -29, 29), (0, 0, 0.2), 31)
        close = camera(
            "Actual DNA texture on enlarged solid cap",
            (13.2, -1.0, 6.9),
            (10.643914, 3.478737, 0.65),
            5.2,
        )
        top = camera("Rigid geographical island plan", (0, 0, 40), (0, 0, 0), 28.5)
    else:
        if cover:
            light(
                "Warm grazing portrait key",
                (-1.6, -2.1, 2.3),
                (0, 0, 0.28),
                160,
                1.8,
                (1.0, 0.87, 0.72),
            )
            light("Soft cool portrait fill", (2.0, -1.4, 1.6), (0, 0, 0.28), 25, 2.6, (0.80, 0.90, 1.0))
            light(
                "Controlled steel rim", (1.3, 1.2, 2.4), (0, 0, 0.35), 110, 1.5, (1.0, 0.95, 0.85)
            )
            light(
                "Silver reflected softbox",
                (-1.15, 0.12, 1.22),
                (0, 0, 0.486),
                7.5,
                1.1,
                (0.92, 0.96, 1.0),
            )
        else:
            light("Large soft grazing key", (-3.5, -3.0, 2.8), (0, -0.7, 0.20), 530, 4.2)
            light("Broad gentle fill", (4.0, -1.2, 4.0), (0, -1.0, 0), 180, 4.0)
            light("Metal reflection strip", (-0.8, 2.2, 3.3), (0, 0, 0.45), 30, 2.0)
        hero = camera("Installation overview", (4.6, -5.5, 5.6), (0.10, -1.2, 0.12), 4.75)
        close = camera(
            "Exact 294-mark solid engraved cap"
            if arguments.cap_treatment == "engraved"
            else "Exact 294-hole cap close view",
            (0.28, -0.40, 1.50),
            (0, 0, 0.43),
            0.62,
        )
        top = camera("Geographical plan view", (0.1365, -1.283, 7.0), (0.1365, -1.283, 0), 3.55)
    if cover:
        portrait_camera = camera(
            "Physical engraved seat portrait", (1.15, -0.12, 1.0), (0, 0, 0.275), 0.82
        )
        portrait_camera.data.type = "PERSP"
        portrait_camera.data.lens = 58
        scene.camera = portrait_camera
        scene.render.resolution_x, scene.render.resolution_y = (
            (640, 800) if arguments.quick else (1600, 2000)
        )
        default_destination = SCENES / (
            "portrait-validation.png" if arguments.quick else "portrait-engraved.png"
        )
        destination = (
            arguments.output.expanduser().resolve() if arguments.output else default_destination
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        scene.render.filepath = str(destination)
        if mesh_signature(meshes) != original_signature:
            raise ValueError("Imported mesh positions, topology or UVs changed during staging.")
        if not arguments.quick:
            bpy.ops.wm.save_as_mainfile(filepath=str(destination.with_suffix(".blend")))
        bpy.ops.render.render(write_still=True)
        if mesh_signature(meshes) != original_signature:
            raise ValueError("Imported mesh positions, topology or UVs changed during rendering.")
        if hashlib.sha256(asset.read_bytes()).hexdigest() != original_asset_sha:
            raise ValueError("Source asset changed during the render.")
        if not arguments.quick:
            provenance = {
                "method": "Actual Blender Cycles CPU render of the imported engraved physical GLB",
                "image_generation_used": False,
                "source_asset": asset.relative_to(ROOT).as_posix(),
                "source_asset_sha256": original_asset_sha,
                "imported_mesh_signature_sha256": original_signature,
                "mesh_positions_topology_and_uvs_unchanged": True,
                "engraved_mark_diameter_mm": 1.5,
                "engraved_substitution_marks": 294,
                "art_direction": "Studio portrait of the generated physical GLB; procedural bark, warm grazing light and restrained silver reflection; no simulated historical architecture.",
                "timber_representation": "Generated source-GLB prototype; not a selected or sourced real timber specimen.",
                "render_only_changes": "Camera, lighting, procedural material shading and studio ground; no scientific geometry changes",
                "camera_position_z_up_m": list(portrait_camera.location),
                "camera_target_z_up_m": [0, 0, 0.275],
                "camera_lens_mm": portrait_camera.data.lens,
                "render_engine": scene.render.engine,
                "samples": scene.cycles.samples,
                "cpu_threads": scene.render.threads,
                "dimensions_px": [scene.render.resolution_x, scene.render.resolution_y],
                "source_is_not_a_site_photograph": True,
            }
            destination.with_suffix(".sources.json").write_text(
                json.dumps(provenance, indent=2) + "\n"
            )
        print(f"Generated physical-seat portrait rendered: {destination}", flush=True)
        return
    IMAGES.mkdir(parents=True, exist_ok=True)
    SCENES.mkdir(parents=True, exist_ok=True)
    scene["source_asset"] = str(asset)
    scene["scientific_geometry"] = (
        "Exact imported GLB; no generated holes or replacement DNA pattern"
    )
    scene["status"] = (
        "Enlarged digital river island edition; separate from physical seat"
        if digital
        else "Illustrative 3 mm prototype; not a structural or venue approval"
    )
    scene.camera = hero
    scene.render.resolution_x, scene.render.resolution_y = 1800, 1500
    if not arguments.quick:
        bpy.ops.wm.save_as_mainfile(filepath=str(SCENES / f"stump-of-life{suffix}.blend"))
    views = (
        [(close, "validation.png", 640, 640)]
        if arguments.quick
        else [
            (hero, f"installation{suffix}.png", 1800, 1200 if digital else 1500),
            (close, f"installation-cap{suffix}.png", 1600, 1600),
            (top, f"installation-top{suffix}.png", 1800, 1000 if digital else 1800),
        ]
    )
    if arguments.view == "cap" and not arguments.quick:
        views = [(close, f"installation-cap{suffix}.png", 1800, 1800)]
    if arguments.view == "top" and not arguments.quick:
        views = [(top, f"installation-top{suffix}.png", 1800, 1000 if digital else 1800)]
    for cam, name, width, height in views:
        scene.camera = cam
        scene.render.resolution_x, scene.render.resolution_y = width, height
        scene.render.filepath = str((SCENES if arguments.quick else IMAGES) / name)
        print(f"Rendering {name} from the exact exported GLB...", flush=True)
        bpy.ops.render.render(write_still=True)
    scene.camera = hero
    scene.render.resolution_x, scene.render.resolution_y = 1800, 1500
    if not arguments.quick:
        bpy.ops.wm.save_as_mainfile(filepath=str(SCENES / f"stump-of-life{suffix}.blend"))
    print("Installation overview, cap close-up and geographical top view complete.", flush=True)


if __name__ == "__main__":
    main()
