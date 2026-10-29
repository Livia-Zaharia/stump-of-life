# Stump of Life contributor guide

- Use `uv`; package sources live in `src/stump_of_life/`.
- Give Python function parameters and returns explicit type annotations; use absolute imports.
- Use Typer for CLI commands. MAFFT is an external dependency, not a bundled binary.
- Preserve every input nucleotide and IUPAC symbol. Never fill unknown DNA from a reference.
- Normalize circles to a common unique anchor before MAFFT. Restore the origin by rotating all
  alignment columns together; never subtract query/reference anchor offsets across indels.
- Distinguish substitutions, insertions, deletions and unknowns. Holes encode called substitutions.
- Full mode retains every alignment column. Variable mode compresses invariant distances and must
  be labelled. Coordinates always name their actual reference; do not silently call them rCRS.
- Ancient samples are observed individuals, not proven ancestors. Study groups are metadata.
- Keep generated runs under gitignored `output/`. Include provenance for public example sequences.
- Run `uv run ruff check .` and `uv run pytest`; install MAFFT for integration tests.
- Retrieval scripts use the locked `research` dependency group (`pysam`, `openpyxl`).
- Keep raw research downloads local under ignored `data/` paths; version manifests and sources.
- Site-map roots encode geographic provenance. Several people at one excavation share one site
  coordinate; never invent separate burial or modern collection coordinates.
- Inspect new SVG output visually. A geometric clearance pass is not structural certification.
