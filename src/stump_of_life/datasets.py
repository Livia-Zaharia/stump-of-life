"""Acquire the fixed public Ko et al. Taiwan mitochondrial sequence panel."""

import hashlib
import json
import re
from datetime import UTC, datetime
from pathlib import Path
from urllib.request import Request, urlopen

from stump_of_life.io import IUPAC, write_fasta

TAIWAN_PANEL_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/KF540505-KF541055?download=false"
TAIWAN_PANEL_DOI = "10.1016/j.ajhg.2014.02.003"
_EXPECTED = tuple(f"KF{number}.1" for number in range(540505, 541056))
_ENA_ID = re.compile(r"ENA\|([A-Z]{1,2}\d+)\|(\1\.[1-9]\d*)")


def parse_ena_fasta(text: str) -> dict[str, str]:
    """Read ENA FASTA into unique versioned accessions and ungapped IUPAC DNA.

    Headers must begin with ``ENA|ACCESSION|ACCESSION.VERSION``. Descriptions
    are discarded; only whitespace and letter case are normalized in sequences.
    This parser does not impose a panel size or a biological sequence length.
    """
    records: dict[str, str] = {}
    identifier: str | None = None
    parts: list[str] = []

    def save() -> None:
        if identifier is None:
            return
        sequence = "".join(parts).upper()
        if not sequence:
            raise ValueError(f"Empty ENA FASTA sequence: {identifier}.")
        invalid = sorted(set(sequence) - IUPAC)
        if invalid:
            raise ValueError(
                f"Invalid ungapped IUPAC DNA in {identifier}: {invalid}. "
                "Sequences are not filled or repaired."
            )
        if identifier in records:
            raise ValueError(f"Duplicate ENA accession: {identifier}.")
        records[identifier] = sequence

    for line_number, raw_line in enumerate(text.splitlines(), 1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            save()
            fields = line[1:].split()
            match = _ENA_ID.fullmatch(fields[0]) if fields else None
            if match is None:
                raise ValueError(
                    f"Invalid ENA FASTA header on line {line_number}; expected "
                    "ENA|ACCESSION|ACCESSION.VERSION with matching accessions."
                )
            identifier = match.group(2)
            parts = []
        elif identifier is None:
            raise ValueError(f"Sequence before an ENA FASTA header on line {line_number}.")
        else:
            parts.append("".join(line.split()))
    save()
    if not records:
        raise ValueError("No ENA FASTA records were returned.")
    return records


def fetch_taiwan_panel(output: Path) -> dict[str, object]:
    """Download, validate and preserve all 551 version-pinned public sequences.

    No genome-length assumption is imposed. Invalid or incomplete responses
    fail before output files are replaced. HTTP and filesystem errors propagate.
    """
    request = Request(
        TAIWAN_PANEL_URL,
        headers={"User-Agent": "stump-of-life (public research data acquisition)"},
    )
    with urlopen(request, timeout=60) as response:
        raw = response.read()
    records = parse_ena_fasta(raw.decode("utf-8"))
    missing = sorted(set(_EXPECTED) - records.keys())
    unexpected = sorted(records.keys() - set(_EXPECTED))
    if missing or unexpected or len(records) != len(_EXPECTED):
        raise ValueError(
            "Expected exactly 551 distinct accessions KF540505.1..KF541055.1; "
            f"received {len(records)}. Missing ({len(missing)}): {missing[:8]}; "
            f"unexpected ({len(unexpected)}): {unexpected[:8]}. "
            "The downloaded accession range is incomplete or has changed."
        )
    ordered = {accession: records[accession] for accession in _EXPECTED}
    raw_sha256 = hashlib.sha256(raw).hexdigest()
    provenance = {
        "dataset": "Ko et al. (2014): Liangdao Man and 550 Taiwanese mitochondrial genomes",
        "source": TAIWAN_PANEL_URL,
        "source_archive": "European Nucleotide Archive",
        "doi": TAIWAN_PANEL_DOI,
        "retrieved_utc": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "count": len(ordered),
        "accession_range": [_EXPECTED[0], _EXPECTED[-1]],
        "raw_file": "ena_original_551.fasta",
        "raw_sha256": raw_sha256,
        "fasta_file": "mitogenomes.fasta",
        "sequences": {
            accession: {
                "length": len(sequence),
                "sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            }
            for accession, sequence in ordered.items()
        },
        "transformations": "FASTA header simplification, whitespace removal, and uppercasing only",
        "sequence_filling": False,
    }
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "ena_original_551.fasta").write_bytes(raw)
    write_fasta(ordered, output / "mitogenomes.fasta")
    (output / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {"count": len(ordered), "output": str(output), "sha256": raw_sha256}
