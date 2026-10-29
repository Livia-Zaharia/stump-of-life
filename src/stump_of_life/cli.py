"""Command-line interface for reproducible DNA-to-pattern runs."""

import json
from enum import StrEnum
from pathlib import Path

import typer

from stump_of_life.alignment import align_circular
from stump_of_life.datasets import fetch_taiwan_panel
from stump_of_life.engraving import export_engraving
from stump_of_life.io import read_fasta
from stump_of_life.metal import export_perforation
from stump_of_life.render import render_alignment_svg

app = typer.Typer(no_args_is_help=True, help="Circular mitochondrial DNA alignment and sculpture patterns.")


class CoordinateMode(StrEnum):
    FULL = "full"
    VARIABLE = "variable"


def _excluded(columns: list[int]) -> list[int]:
    if any(column < 1 for column in columns):
        raise typer.BadParameter("Excluded alignment columns are one-based positive numbers")
    return [column - 1 for column in columns]


@app.command("fetch-taiwan")
def fetch_taiwan(
    output: Path = typer.Option(Path("data/modern/taiwan-panel")),
) -> None:
    """Download the public Ko 2014 panel: one ancient and 550 modern mitogenomes."""
    typer.echo(json.dumps(fetch_taiwan_panel(output), indent=2))


@app.command()
def align(
    input: Path = typer.Argument(..., exists=True, dir_okay=False),
    reference_id: str = typer.Option(..., help="Exact FASTA ID; its original base 1 defines coordinates."),
    output: Path = typer.Option(Path("output/alignment")),
    mafft: str = typer.Option("mafft", help="MAFFT executable or path."),
    mafft_binaries: str | None = typer.Option(None, help="Relocated MAFFT library directory."),
    threads: int = typer.Option(2, min=1),
    anchor_size: int = typer.Option(31, min=3),
) -> None:
    """Normalize complete DNA circles, align with MAFFT, and restore reference coordinates."""
    summary = align_circular(
        read_fasta(input), reference_id, output, mafft, threads, mafft_binaries, anchor_size,
    )
    typer.echo(json.dumps(summary, indent=2))
    typer.echo(f"Alignment and coordinate tables: {output.resolve()}")


@app.command()
def render(
    input: Path = typer.Argument(..., exists=True, dir_okay=False),
    reference_id: str = typer.Option(...),
    output: Path = typer.Option(Path("output/alignment.svg")),
    mode: CoordinateMode = typer.Option(CoordinateMode.FULL),
    diameter: float = typer.Option(500, min=1, help="SVG diameter in millimetres."),
    select: list[str] | None = typer.Option(None, help="Outer sample ID; repeat to select several."),
    exclude_column: list[int] | None = typer.Option(None, help="One-based column to suppress; repeatable."),
) -> None:
    """Draw an aligned circular comparison, preserving unknown and indel states."""
    summary = render_alignment_svg(
        read_fasta(input, aligned=True), reference_id, output,
        selected_ids=select, coordinate_mode=mode.value, diameter_mm=diameter,
        excluded_columns=_excluded(exclude_column or []),
    )
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    typer.echo(f"Circular alignment: {output.resolve()}")


@app.command()
def engrave(
    input: Path = typer.Argument(..., exists=True, dir_okay=False),
    reference_id: str = typer.Option(...),
    output: Path = typer.Option(Path("output/engraving")),
    mode: CoordinateMode = typer.Option(CoordinateMode.VARIABLE),
    diameter: float = typer.Option(500, min=1),
    mark_diameter: float = typer.Option(1.2, min=0.01),
    minimum_gap: float = typer.Option(0.2, min=0),
    select: list[str] | None = typer.Option(None),
    exclude_column: list[int] | None = typer.Option(None),
) -> None:
    """Export a solid cap with surface marks; CUT contains only the outer edge."""
    summary = export_engraving(
        read_fasta(input, aligned=True), reference_id, output,
        diameter_mm=diameter, mark_diameter_mm=mark_diameter, min_gap_mm=minimum_gap,
        coordinate_mode=mode.value, selected_ids=select,
        excluded_columns=_excluded(exclude_column or []),
    )
    typer.echo(json.dumps(summary, indent=2))
    if not summary["passes_requested_clearance"]:
        raise typer.Exit(code=2)


@app.command()
def perforate(
    input: Path = typer.Argument(..., exists=True, dir_okay=False),
    reference_id: str = typer.Option(...),
    output: Path = typer.Option(Path("output/metal")),
    mode: CoordinateMode = typer.Option(CoordinateMode.VARIABLE),
    diameter: float = typer.Option(500, min=1),
    hole_diameter: float = typer.Option(5, min=0.01),
    minimum_web: float = typer.Option(1.5, min=0),
    select: list[str] | None = typer.Option(None),
    exclude_column: list[int] | None = typer.Option(None),
) -> None:
    """Export millimetre cut/engrave SVG and check hole and edge clearances."""
    summary = export_perforation(
        read_fasta(input, aligned=True), reference_id, output,
        diameter_mm=diameter, hole_diameter_mm=hole_diameter,
        min_web_mm=minimum_web, coordinate_mode=mode.value,
        selected_ids=select, excluded_columns=_excluded(exclude_column or []),
    )
    typer.echo(json.dumps(summary, indent=2))
    typer.echo(f"Design geometry: {output.resolve()}")
    if not summary["passes_requested_clearance"]:
        typer.echo("Requested clearances FAIL. Review geometry.json before fabrication.", err=True)
        raise typer.Exit(code=2)
