"""Fetch the two small public Suogang BAMs and retain their observed MT reads.

Run with ``uv run --with pysam python scripts/fetch_suogang.py``.
This does not call a mitochondrial consensus or replace missing data with a reference.
"""

from __future__ import annotations

import csv
import hashlib
import html
import json
import re
import time
import urllib.error
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pysam
import typer

ARCHIVE = "https://download.cncb.ac.cn/gsa-human/HRA000123"
PAPER = "https://doi.org/10.1126/science.aba0909"
SUPPLEMENT = (
    "https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/"
    "YangSciecne2020Supplement.pdf"
)
RUNS = [
    {
        "individual_id": "SuogangB1",
        "run_accession": "HRR051945",
        "sample_accession": "HRS029483",
        "individual_accession": "HRI029483",
        "experiment_accession": "HRX035744",
        "expected_bytes": 3651852,
        "expected_md5": "cc255bd293befc1aa8e05b4522e16a00",
    },
    {
        "individual_id": "SuogangB3",
        "run_accession": "HRR051946",
        "sample_accession": "HRS029484",
        "individual_accession": "HRI029484",
        "experiment_accession": "HRX035745",
        "expected_bytes": 2212862,
        "expected_md5": "f2b9985c694cc182c95d0dfc49ee01f6",
    },
]


def digest(path: Path, algorithm: str = "sha256") -> str:
    """Hash a local file without loading it into memory."""
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def fetch(url: str, destination: Path, expected_bytes: int | None = None) -> None:
    """Download atomically, retrying transient connection failures on these small files."""
    if destination.exists() and (
        expected_bytes is None or destination.stat().st_size == expected_bytes
    ):
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "Stump-of-Life/0.1"})
            with urllib.request.urlopen(request, timeout=40) as response:
                with temporary.open("wb") as handle:
                    while block := response.read(1024 * 1024):
                        handle.write(block)
            if expected_bytes is not None and temporary.stat().st_size != expected_bytes:
                raise OSError(f"Unexpected size for {url}: {temporary.stat().st_size}")
            temporary.replace(destination)
            return
        except (OSError, urllib.error.URLError):
            if attempt == 3:
                raise
            time.sleep(attempt + 1)


def inspect_and_extract(source: Path, output: Path) -> dict[str, Any]:
    """Extract deposited MT alignments unchanged and count observable base coverage."""
    with pysam.AlignmentFile(source, "rb") as bam:
        reference_names = bam.references
        mitochondrial = [name for name in reference_names if name in {"MT", "chrM", "M"}]
        if len(mitochondrial) != 1:
            raise ValueError(f"Expected one mitochondrial contig in {source}")
        contig = mitochondrial[0]
        length = bam.get_reference_length(contig)
        depth_all = [0] * length
        depth_q30 = [0] * length
        read_count = 0
        read_count_q30 = 0
        alignments_by_contig: Counter[str] = Counter()
        library_names: set[str] = set()
        mt_header = dict(bam.header.to_dict())
        references = mt_header["SQ"]
        with pysam.AlignmentFile(output, "wb", header=bam.header) as target:
            for read in bam.fetch(until_eof=True):
                alignments_by_contig[read.reference_name or "unmapped"] += 1
                if read.is_unmapped or read.reference_name != contig:
                    continue
                target.write(read)
                read_count += 1
                if read.mapping_quality >= 30:
                    read_count_q30 += 1
                if read.has_tag("RG"):
                    library_names.add(read.get_tag("RG"))
                qualities = read.query_qualities
                for query_position, reference_position in read.get_aligned_pairs(matches_only=True):
                    if query_position is None or reference_position is None:
                        continue
                    depth_all[reference_position] += 1
                    if (
                        read.mapping_quality >= 30
                        and qualities is not None
                        and qualities[query_position] >= 30
                        and not read.is_secondary
                        and not read.is_supplementary
                        and not read.is_duplicate
                        and not read.is_qcfail
                    ):
                        depth_q30[reference_position] += 1
        pysam.index(str(output))
    coverage_path = output.with_suffix(".coverage.csv")
    with coverage_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["reference_position", "all_aligned_bases", "mapq30_baseq30_primary_bases"])
        writer.writerows((i + 1, depth_all[i], depth_q30[i]) for i in range(length))

    def summarize(depth: list[int]) -> dict[str, float | int]:
        return {
            "mean_depth": sum(depth) / length,
            **{f"bases_at_least_{n}x": sum(value >= n for value in depth) for n in (1, 3, 5)},
        }

    return {
        "source_reference_name": contig,
        "source_reference_length": length,
        "source_reference_record": next(row for row in references if row["SN"] == contig),
        "reference_identity": (
            "BAM @PG states hg19_1000g; MT length alone does not prove rCRS identity. "
            "No sequence reference is supplied in this deposit or substituted here."
        ),
        "mt_aligned_reads": read_count,
        "total_deposited_alignment_records": sum(alignments_by_contig.values()),
        "deposited_alignment_records_by_contig": dict(alignments_by_contig),
        "mt_reads_mapping_quality_at_least_30": read_count_q30,
        "all_aligned_base_coverage": summarize(depth_all),
        "mapq30_baseq30_primary_coverage": summarize(depth_q30),
        "read_group_tags_observed": sorted(library_names),
        "mt_bam": str(output),
        "mt_bam_sha256": digest(output),
        "coverage_csv": str(coverage_path),
        "consensus_generated": False,
        "mt_sequence_status": (
            "No mitochondrial reads in deposited BAM" if read_count == 0 else "Observed reads only"
        ),
        "damage_filter_applied_by_this_script": False,
    }


def main(
    output: Path = Path("data/ancient/suogang"),
    max_download_mb: float = 25,
) -> None:
    """Download just two Suogang individuals; never scan or bulk-fetch the study's BAMs."""
    if sum(int(row["expected_bytes"]) for row in RUNS) > max_download_mb * 1_000_000:
        raise typer.BadParameter("The two BAMs exceed the requested download cap")
    output.mkdir(parents=True, exist_ok=True)
    checksum_file = output / "archive_md5sum.txt"
    fetch(f"{ARCHIVE}/md5sum.txt", checksum_file)
    checksums = checksum_file.read_text()
    records: list[dict[str, Any]] = []
    for sample in RUNS:
        row = dict(sample)
        run = str(row["run_accession"])
        metadata_url = f"https://ngdc.cncb.ac.cn/gsa-human/browse/runDetail/{run}"
        metadata_path = output / f"{run}.html"
        fetch(metadata_url, metadata_path)
        title = re.search(
            r"<th>\s*title\s*</th>\s*<td>(.*?)</td>",
            metadata_path.read_text(),
            re.DOTALL,
        )
        if title is None or html.unescape(title.group(1)).strip() != row["individual_id"]:
            raise ValueError(f"Official run page does not identify expected sample: {run}")
        for label, field in (
            ("Individual", "individual_accession"),
            ("Sample", "sample_accession"),
            ("Experiment", "experiment_accession"),
        ):
            match = re.search(
                rf"<th>\s*{label}\s*</th>\s*<td>(.*?)</td>",
                metadata_path.read_text(),
                re.DOTALL,
            )
            value = re.sub(r"<[^>]+>", "", match.group(1)).strip() if match else None
            if value != row[field]:
                raise ValueError(f"Unexpected {label} accession on official page: {run}")
        bam_url = f"{ARCHIVE}/{run}/{run}.bam"
        source = output / f"{run}.bam"
        fetch(bam_url, source, int(row["expected_bytes"]))
        md5 = digest(source, "md5")
        if md5 != row["expected_md5"] or md5 not in checksums:
            raise ValueError(f"Archive MD5 mismatch: {run}")
        row.update(
            {
                "source_url": bam_url,
                "metadata_url": metadata_url,
                "metadata_sha256": digest(metadata_path),
                "bytes": source.stat().st_size,
                "md5_verified": True,
                "sha256": digest(source),
                "site": "Suogang, Penghu archipelago",
                "latitude": 23.5,
                "longitude": 119.6,
                "location_precision": "Rounded study Table 1 coordinates, not burial coordinates",
                "date_reported_cal_bp": "4800–4300",
                "date_method": "Contextual associated-shell radiocarbon dates; bone not directly dated",
                "site_associated_shell_ranges_cal_bp": [[4633, 4287], [4793, 4407]],
                "individual_to_shell_range_pairing": "Not independently verified; retain both site ranges",
                "contamination": (
                    "Yang2020 SI reports contamination; published population analyses used "
                    "damage-restricted fragments. These downloaded BAMs are not asserted clean."
                ),
            }
        )
        row.update(inspect_and_extract(source, output / f"{row['individual_id']}.MT.bam"))
        records.append(row)
        typer.echo(f"{row['individual_id']}: {row['mt_aligned_reads']} deposited MT reads")
    manifest = {
        "study_accession": "HRA000123",
        "bioproject": "PRJCA002264",
        "publication_doi": "10.1126/science.aba0909",
        "paper_url": PAPER,
        "supplement_url": SUPPLEMENT,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "acquisition": "Two small public mixed-reference BAMs, archive MD5 verified; MT extracted unchanged",
        "total_bam_download_bytes": sum(int(row["bytes"]) for row in records),
        "complete_mitogenome_fasta_identified": False,
        "consensus_generated": False,
        "limitations": [
            "Coverage measures deposited alignments; it is not an authenticity or contamination estimate.",
            "No ancient individual is assumed to be an ancestor of another sample or modern group.",
            "The displayed dates are archaeological contexts, not exact individual death dates.",
            "No missing bases are filled from a reference.",
        ],
        "individuals": records,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (output / "individuals.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "individual_id", "run_accession", "sample_accession", "site", "latitude",
                "longitude", "date_reported_cal_bp", "date_method", "mt_aligned_reads",
                "mt_sequence_status", "source_url", "metadata_url",
            ],
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    typer.run(main)
