"""Fetch published Liangdao-1 mtDNA and mitochondrial-only Liangdao-2 reads.

Run from the repository root with ``uv run --with pysam python
scripts/fetch_liangdao.py``. Coordinate-sorted ENA BAMs support byte ranges but
have no public BAI sidecar. Bounded BGZF probes locate the mitochondrial band;
only that band is downloaded. No consensus is invented from the reads.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import struct
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.request import Request, urlopen

import pysam

from stump_of_life.datasets import parse_ena_fasta
from stump_of_life.io import write_fasta

OUTPUT = Path(__file__).resolve().parents[1] / "data" / "ancient" / "liangdao"
PUBLICATION_DOI = "10.1073/pnas.2026132118"
PUBLICATION_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC8020671/"
REPORT_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB43078"
    "&result=read_run&fields=study_accession,sample_accession,sample_alias,run_accession,"
    "submitted_ftp,submitted_bytes,submitted_md5,submitted_format&format=json"
)
FASTA_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/KF540505.1?download=false"
MAGIC = b"\x1f\x8b\x08\x04"
PROBE_BYTES = 262144
MAX_SLICE_BYTES = 32 * 1024 * 1024
USER_AGENT = "stump-of-life (bounded public mitochondrial data acquisition)"


def download(url: str) -> bytes:
    with urlopen(Request(url, headers={"User-Agent": USER_AGENT}), timeout=60) as response:
        return response.read()


def ranged(url: str, start: int, end: int, audit: list[dict[str, object]]) -> bytes:
    """Require range support; never fall back to a whole-genome download."""
    if end - start + 1 > MAX_SLICE_BYTES:
        raise ValueError("Refusing a mitochondrial source slice larger than 32 MiB.")
    request = Request(url, headers={"User-Agent": USER_AGENT, "Range": f"bytes={start}-{end}"})
    with urlopen(request, timeout=60) as response:
        if response.status != 206:
            raise ValueError("Server did not return HTTP 206; refusing the full BAM.")
        content_range = response.headers.get("Content-Range", "")
        if not content_range.startswith(f"bytes {start}-{end}/"):
            raise ValueError(f"Unexpected Content-Range: {content_range!r}.")
        raw = response.read(end - start + 2)
    if len(raw) != end - start + 1:
        raise ValueError("Incomplete or oversized HTTP range response.")
    audit.append({"start": start, "end": end, "bytes": len(raw), "sha256": sha256(raw)})
    return raw


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def decode_blocks(raw: bytes) -> tuple[int, bytes]:
    """Decode complete BGZF members after the first boundary in a byte range."""
    first = raw.find(MAGIC)
    if first < 0:
        raise ValueError("No BGZF block boundary in the range probe.")
    offset = first
    decoded: list[bytes] = []
    while offset + 18 <= len(raw):
        if raw[offset : offset + 4] != MAGIC or raw[offset + 12 : offset + 14] != b"BC":
            raise ValueError("Invalid BGZF member at a detected block boundary.")
        size = struct.unpack_from("<H", raw, offset + 16)[0] + 1
        if offset + size > len(raw):
            break
        decoded.append(gzip.decompress(raw[offset : offset + size]))
        offset += size
    if not decoded:
        raise ValueError("No complete BGZF members in the range probe.")
    return first, b"".join(decoded)


def parse_header(decoded: bytes) -> tuple[bytes, list[tuple[str, int]]]:
    if decoded[:4] != b"BAM\x01":
        raise ValueError("The source is not BAM.")
    text_length = struct.unpack_from("<i", decoded, 4)[0]
    text = decoded[8 : 8 + text_length].decode("ascii")
    if "SO:coordinate" not in text.splitlines()[0]:
        raise ValueError("Mitochondrial range extraction requires a coordinate-sorted BAM.")
    cursor = 8 + text_length
    reference_count = struct.unpack_from("<i", decoded, cursor)[0]
    cursor += 4
    references: list[tuple[str, int]] = []
    for _ in range(reference_count):
        name_length = struct.unpack_from("<i", decoded, cursor)[0]
        cursor += 4
        name = decoded[cursor : cursor + name_length - 1].decode("ascii")
        cursor += name_length
        length = struct.unpack_from("<i", decoded, cursor)[0]
        cursor += 4
        references.append((name, length))
    return decoded[:cursor], references


def record_info(
    data: bytes, offset: int, references: list[tuple[str, int]]
) -> tuple[int, int, int] | None:
    """Validate a possible BAM record boundary before accepting a range probe."""
    if offset + 36 > len(data):
        return None
    size, reference_id, position, packed, flags_cigar, sequence_length = struct.unpack_from(
        "<iiiiii", data, offset
    )
    end = offset + 4 + size
    if not 32 <= size <= 1000000 or end > len(data):
        return None
    if not -1 <= reference_id < len(references):
        return None
    if reference_id >= 0 and not -1 <= position < references[reference_id][1]:
        return None
    name_length = packed & 255
    cigar_count = flags_cigar & 65535
    if not 2 <= name_length <= 255 or not 0 <= sequence_length <= 100000:
        return None
    required = 32 + name_length + 4 * cigar_count + (sequence_length + 1) // 2 + sequence_length
    if required > size:
        return None
    name = data[offset + 36 : offset + 36 + name_length]
    if (
        len(name) != name_length
        or name[-1] != 0
        or not all(32 <= value <= 126 for value in name[:-1])
    ):
        return None
    return end, reference_id, position


def record_start(data: bytes, references: list[tuple[str, int]]) -> int:
    """Find a boundary by requiring three consecutive structurally valid records."""
    for offset in range(min(65536, len(data))):
        cursor = offset
        for _ in range(3):
            info = record_info(data, cursor, references)
            if info is None:
                break
            cursor = info[0]
        else:
            return offset
    raise ValueError("Could not establish a reliable BAM record boundary in the probe.")


def coverage_summary(reads: list[Any], reference_length: int) -> dict[str, object]:
    """Describe observed support without calling a mitochondrial consensus."""
    raw_depth = [0] * reference_length
    quality_depth = [0] * reference_length
    for read in reads:
        eligible = read.mapping_quality >= 30 and not (
            read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail
        )
        qualities = read.query_qualities
        for query, reference in read.get_aligned_pairs(matches_only=True):
            if not 0 <= reference < reference_length:
                continue
            raw_depth[reference] += 1
            if eligible and qualities is not None and qualities[query] >= 30:
                quality_depth[reference] += 1
    return {
        "reference_length": reference_length,
        "raw_mean_depth": sum(raw_depth) / reference_length,
        "raw_covered_positions": sum(depth > 0 for depth in raw_depth),
        "quality_filter": "MQ >= 30, BQ >= 30; exclude secondary, supplementary, duplicate, QC-failed reads",
        "quality_filtered_mean_depth": sum(quality_depth) / reference_length,
        "quality_filtered_covered_positions": sum(depth > 0 for depth in quality_depth),
        "quality_filtered_zero_depth_positions": sum(depth == 0 for depth in quality_depth),
        "quality_filtered_depth_at_least_three_positions": sum(
            depth >= 3 for depth in quality_depth
        ),
        "read_end_trimming": "none; this is a descriptive coverage audit, not damage-aware consensus calling",
        "consensus_generated": False,
    }


def extract_mitochondrial_reads(row: dict[str, str], folder: Path) -> dict[str, object]:
    url = "https://" + row["submitted_ftp"]
    size = int(row["submitted_bytes"])
    audit: list[dict[str, object]] = []
    head = ranged(url, 0, PROBE_BYTES - 1, audit)
    _, decoded = decode_blocks(head)
    header, references = parse_header(decoded)
    mitochondrial = [index for index, (name, _) in enumerate(references) if name in {"MT", "chrM"}]
    if len(mitochondrial) != 1:
        raise ValueError("Expected exactly one mitochondrial reference in BAM header.")
    target = mitochondrial[0]

    def probe(offset: int) -> int:
        raw = ranged(url, offset, min(size - 1, offset + PROBE_BYTES - 1), audit)
        _, data = decode_blocks(raw)
        info = record_info(data, record_start(data, references), references)
        assert info is not None
        return info[1] if info[1] >= 0 else len(references)

    def boundary(reference_id: int) -> tuple[int, int]:
        low, high = 0, size - PROBE_BYTES
        while high - low > PROBE_BYTES:
            middle = (low + high) // 2
            if probe(middle) >= reference_id:
                high = middle
            else:
                low = middle
        return low, high

    before, _ = boundary(target)
    _, after = boundary(target + 1)
    start = max(0, before - PROBE_BYTES)
    end = min(size - 1, after + 2 * PROBE_BYTES)
    raw = ranged(url, start, end, audit)
    _, data = decode_blocks(raw)
    cursor = record_start(data, references)
    mitochondrial_records: list[bytes] = []
    preceding = following = False
    previous_coordinate = (-1, -1)
    while (info := record_info(data, cursor, references)) is not None:
        next_cursor, reference_id, position = info
        sort_id = reference_id if reference_id >= 0 else len(references)
        if (sort_id, position) < previous_coordinate:
            raise ValueError("Source coordinates are not monotonically sorted.")
        previous_coordinate = (sort_id, position)
        preceding |= 0 <= reference_id < target
        following |= reference_id > target or reference_id == -1
        if reference_id == target:
            mitochondrial_records.append(data[cursor:next_cursor])
        cursor = next_cursor
    if not mitochondrial_records or not preceding or not following:
        raise ValueError(
            "Range slice does not bracket all mitochondrial reads; refusing incomplete extraction."
        )
    run = row["run_accession"]
    bam_path = folder / f"{run}.MT.bam"
    with pysam.BGZFile(str(bam_path), "wb") as bam:
        bam.write(header)
        for record in mitochondrial_records:
            bam.write(record)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        checked = list(bam.fetch(until_eof=True))
        if len(checked) != len(mitochondrial_records) or any(
            read.reference_id != target for read in checked
        ):
            raise ValueError("Repacked mitochondrial BAM failed validation.")
    pysam.index(str(bam_path))
    source_slice = folder / f"{run}.MT-source-range.bgzf-part"
    source_slice.write_bytes(raw)
    (folder / f"{run}.head.bgzf-part").write_bytes(head)
    return {
        "run_accession": run,
        "sample_accession": row["sample_accession"],
        "sample_alias": row["sample_alias"],
        "source_url": url,
        "source_whole_file_bytes": size,
        "source_whole_file_md5_declared_by_ena": row["submitted_md5"],
        "whole_file_downloaded": False,
        "source_slice_start": start,
        "source_slice_end": end,
        "source_slice_file": source_slice.name,
        "source_slice_sha256": sha256(raw),
        "downloaded_bytes_including_probes": sum(int(item["bytes"]) for item in audit),
        "range_audit": audit,
        "mitochondrial_reference": references[target][0],
        "mitochondrial_reference_length": references[target][1],
        "read_count": len(checked),
        "coverage_summary": coverage_summary(checked, references[target][1]),
        "bam_file": bam_path.name,
        "bam_sha256": sha256(bam_path.read_bytes()),
        "consensus_generated": False,
        "transformation": "Original BAM header and all MT alignment records repacked unchanged into BGZF BAM",
    }


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    report_raw = download(REPORT_URL)
    (OUTPUT / "ena_run_report.json").write_bytes(report_raw)
    runs = json.loads(report_raw)
    raw_fasta = download(FASTA_URL)
    published = parse_ena_fasta(raw_fasta.decode("utf-8"))
    if set(published) != {"KF540505.1"}:
        raise ValueError("Unexpected published Liangdao-1 accession response.")
    (OUTPUT / "Liangdao1_KF540505.1.ena-original.fasta").write_bytes(raw_fasta)
    write_fasta(published, OUTPUT / "Liangdao1_KF540505.1.fasta")
    for accession in ("SAMEA8393779", "SAMEA8393780"):
        (OUTPUT / f"{accession}.xml").write_bytes(
            download(f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}")
        )
    selected = next(row for row in runs if row["run_accession"] == "ERR5525699")
    print("Locating Liangdao-2 mitochondrial reads with bounded HTTP ranges...", flush=True)
    mito = extract_mitochondrial_reads(selected, OUTPUT)
    sequence = published["KF540505.1"]
    manifest = {
        "dataset": "Liangdao archaeological individuals, Matsu archipelago",
        "retrieved_utc": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "study_accession": "PRJEB43078",
        "study_report_url": REPORT_URL,
        "study_report_sha256": sha256(report_raw),
        "publication_doi": PUBLICATION_DOI,
        "publication_url": PUBLICATION_URL,
        "site_coordinates": {
            "latitude": 26.344958,
            "longitude": 120.227029,
            "role": "Official archaeological heritage-site location; not individual burial GPS",
            "source": "https://tcmb.culture.tw/zh-tw/detail?id=20210618000002&indexCode=BOCH_CountryCulture_21",
            "published_study_latitude": 26.3,
            "published_study_longitude": 120.2,
        },
        "individuals": [
            {
                "individual": "Liangdao-1",
                "run_alias": "tai001",
                "sex_reported": "male",
                "date_cal_bp_min": 8060,
                "date_cal_bp_max": 8320,
                "published_mitogenome_accession": "KF540505.1",
                "published_mitogenome_source": FASTA_URL,
                "published_mitogenome_doi": "10.1016/j.ajhg.2014.02.003",
                "published_mitogenome_length": len(sequence),
                "published_mitogenome_sequence_sha256": sha256(sequence.encode("ascii")),
                "published_mitogenome_raw_sha256": sha256(raw_fasta),
                "complete_published_sequence_available": True,
                "later_sequencing_run": "ERR5525679",
            },
            {
                "individual": "Liangdao-2",
                "run_alias": "tai002",
                "sex_reported": "female",
                "date_cal_bp_min": 7510,
                "date_cal_bp_max": 7590,
                "published_haplogroup": "R9 (study's reported classification)",
                "complete_published_sequence_available": False,
                "status": "Mitochondrial reads extracted; no newly inferred consensus or missing-base filling",
                "mitochondrial_reads": mito,
            },
        ],
        "date_convention": "calibrated years before present (BP), where present is 1950 CE",
        "identity_source": "Larena et al. 2021 SI sections 7.2-7.4: tai001=Liangdao-1; tai002=Liangdao-2",
        "metadata_note": "ENA sample tissue_lib says saliva; ancient identity and dates follow the archaeological publication instead of this inconsistent deposition field.",
        "geography_note": "Liangdao is within the Matsu archipelago; Matsu is not a separate archaeological individual or cohort.",
        "interpretation": "Independent archaeological individuals; no direct maternal ancestry claim",
        "consensus_note": "Do not treat raw reads as a validated complete mitogenome. Contamination-aware consensus calling and explicit masks are a separate analysis.",
    }
    (OUTPUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "published_sequences": 1,
                "Liangdao2_mt_reads": mito["read_count"],
                "downloaded_bytes_including_probes": mito["downloaded_bytes_including_probes"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
