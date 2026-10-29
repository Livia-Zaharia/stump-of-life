"""Fetch Wang 2021 site metadata and indexed mitochondrial read subsets.

Run: uv run --with openpyxl --with pysam python scripts/fetch_hanben_gongguan.py
Source nuclear BAMs are never retained. An unindexed source is streamed through
a mitochondrial filter once; otherwise only indexed ranges are read. No
consensus is inferred.
"""

import csv
import hashlib
import io
import json
import re
import subprocess
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, build_opener

import openpyxl
import pysam
import typer

STUDY = "PRJEB42781"
PAPER = "https://pmc.ncbi.nlm.nih.gov/articles/PMC7993749/"
TABLE_URL = (
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/7993749/bin/"
    "NIHMS1671260-supplement-1671260_OnlineTable1.xlsx"
)
FIELDS = (
    "run_accession,sample_accession,secondary_sample_accession,sample_alias,sample_title,"
    "submitted_ftp,submitted_md5,submitted_bytes,fastq_ftp,fastq_bytes"
)
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urlencode({
    "accession": STUDY, "result": "read_run", "fields": FIELDS, "format": "tsv",
})
DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "data/ancient/hanben_gongguan"


def digest(path: Path, algorithm: str = "sha256") -> str:
    checksum = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def fetch(url: str) -> bytes:
    """Retry transient public-download failures without suppressing final errors."""
    for attempt in range(3):
        try:
            with build_opener().open(url, timeout=120) as response:
                return response.read()
        except OSError:
            if attempt == 2:
                raise
            time.sleep(attempt + 1)
    raise RuntimeError("Unreachable download state")


def fetch_table(path: Path) -> None:
    if path.exists() and path.read_bytes().startswith(b"PK"):
        return
    content = fetch(TABLE_URL)
    if not content.startswith(b"PK"):
        # Perform the public viewer's documented client proof-of-work, rather
        # than treating its HTML preparation page as a downloaded spreadsheet.
        text = content.decode("utf-8")
        challenge_match = re.search(r'POW_CHALLENGE = "([^"]+)"', text)
        difficulty_match = re.search(r'POW_DIFFICULTY = "(\d+)"', text)
        if challenge_match is None or difficulty_match is None:
            raise ValueError("PMC did not return the expected public XLSX or download preparation page")
        challenge = challenge_match[1]
        difficulty = int(difficulty_match[1])
        if not 1 <= difficulty <= 5:
            raise ValueError("Unexpected PMC preparation-work difficulty")
        nonce = 0
        while not hashlib.sha256(f"{challenge}{nonce}".encode()).hexdigest().startswith(
            "0" * difficulty,
        ):
            nonce += 1
        request = Request(TABLE_URL, headers={
            "Cookie": f"cloudpmc-viewer-pow={challenge},{nonce}",
        })
        with build_opener().open(request, timeout=120) as response:
            content = response.read()
    if not content.startswith(b"PK"):
        raise ValueError("PMC response is not an XLSX file; refusing to save HTML as metadata")
    path.write_bytes(content)


def individual_metadata(table: Path) -> list[dict[str, Any]]:
    workbook = openpyxl.load_workbook(table, read_only=True, data_only=True)
    individuals: list[dict[str, Any]] = []
    for number, row in enumerate(workbook.active.values, 1):
        if len(row) <= 35 or row[13] != "Taiwan":
            continue
        if not any(site in str(row[11]) for site in ("Hanben", "Gongguan")):
            continue
        individuals.append({
            "individual_id": row[0], "site": row[11], "skeletal_code": row[1],
            "date_method": row[4], "date_mean_bp": row[5], "date_sd_bp": row[6],
            "date_reported": row[7], "group_id": row[9],
            "included_in_study_analysis": row[10], "latitude": row[14], "longitude": row[15],
            "data_type": row[16], "library_count": row[17], "sex": row[20],
            "reported_mt_coverage": row[25], "reported_mt_haplogroup": row[26],
            "reported_mt_match_to_consensus": row[27], "library_type": row[32],
            "library_ids": row[33], "quality_assessment": row[35], "source_table_row": number,
        })
    workbook.close()
    if len(individuals) != len({row["individual_id"] for row in individuals}):
        raise ValueError("The primary table has repeated individual identifiers; mapping needs review")
    return individuals


def passing_read(read: pysam.AlignedSegment) -> bool:
    return read.mapping_quality >= 30 and not (read.flag & 0xF04)


def stream_unindexed_mt(bam_url: str, target: Path, expected_md5: str, expected_bytes: int) -> None:
    """Filter a public unindexed BAM while hashing its stream, without saving it."""
    expression = 'rname == "MT" || rname == "chrM" || rname == "M" || rname == "chrMT"'
    checksum = hashlib.md5()
    received = 0
    with subprocess.Popen(
        ["samtools", "view", "-b", "-o", str(target), "-e", expression, "-"],
        stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
    ) as process:
        assert process.stdin is not None
        assert process.stderr is not None
        try:
            with build_opener().open(bam_url, timeout=120) as response:
                while block := response.read(1024 * 1024):
                    checksum.update(block)
                    received += len(block)
                    process.stdin.write(block)
            process.stdin.close()
            diagnostic = process.stderr.read().decode("utf-8", errors="replace")
            return_code = process.wait()
        except BaseException:
            # Terminate the local consumer if the public stream fails; no
            # partial subset is promoted to a successful acquisition.
            process.kill()
            process.wait()
            target.unlink(missing_ok=True)
            raise
    if return_code:
        target.unlink(missing_ok=True)
        raise RuntimeError(f"samtools mitochondrial stream filter failed: {diagnostic}")
    if received != expected_bytes or checksum.hexdigest() != expected_md5:
        target.unlink(missing_ok=True)
        raise ValueError("The complete streamed source BAM failed ENA size/MD5 verification")


def extract_run(
    individual: dict[str, Any], run: dict[str, str], output: Path,
) -> dict[str, Any]:
    files = list(zip(
        run["submitted_ftp"].split(";"), run["submitted_md5"].split(";"),
        run["submitted_bytes"].split(";"), strict=True,
    ))
    bams = [item for item in files if item[0].endswith(".bam")]
    if len(bams) != 1:
        raise ValueError(f"Expected one submitted BAM for {run['run_accession']}; found {len(bams)}")
    bam_ftp, bam_md5, bam_bytes = bams[0]
    indices = [item for item in files if item[0] in {bam_ftp + ".bai", bam_ftp[:-4] + ".bai"}]
    if len(indices) > 1:
        raise ValueError(f"Multiple public BAI candidates need review: {run['run_accession']}")
    bam_url = "https://" + bam_ftp
    stem = f"{individual['individual_id']}_{run['run_accession']}"
    index_path = output / "source_indexes" / f"{stem}.source.bam.bai"
    if indices:
        index_ftp, index_md5, index_bytes = indices[0]
        index_url = "https://" + index_ftp
        if not index_path.exists() or digest(index_path, "md5") != index_md5:
            content = fetch(index_url)
            if len(content) != int(index_bytes) or hashlib.md5(content).hexdigest() != index_md5:
                raise ValueError(f"Downloaded source index failed ENA size/MD5 verification: {stem}")
            index_path.write_bytes(content)
        acquisition = {
            "acquisition_method": "indexed_remote_ranges",
            "source_bam_fully_downloaded": False, "source_bam_md5_verified": False,
            "source_bam_retained": False, "source_index_url": index_url,
            "source_index_bytes": int(index_bytes), "source_index_md5": index_md5,
            "source_index_md5_verified": True, "source_index_sha256": digest(index_path),
            "source_index_path": str(index_path),
        }
    else:
        acquisition = {
            "acquisition_method": "complete_stream_mt_filter_no_source_bam_retained",
            "source_bam_fully_downloaded": True, "source_bam_md5_verified": True,
            "source_bam_retained": False, "source_index_url": None,
            "source_index_bytes": None, "source_index_md5": None,
            "source_index_md5_verified": None, "source_index_sha256": None,
            "source_index_path": None,
        }
    subset = output / "mt_reads" / f"{stem}.mt.bam"
    pending = subset.with_suffix(".bam.partial")
    acquisition_path = output / "mt_reads" / f"{stem}.acquisition.json"
    if not subset.exists() or (not indices and not acquisition_path.exists()):
        if indices:
            with pysam.AlignmentFile(bam_url, "rb", index_filename=str(index_path)) as source:
                candidates = [name for name in source.references
                              if name in {"MT", "chrM", "M", "chrMT"}]
                if len(candidates) != 1:
                    raise ValueError(f"Mitochondrial reference is not uniquely identifiable: {stem}")
                contig = candidates[0]
                with pysam.AlignmentFile(str(pending), "wb", template=source) as target:
                    for read in source.fetch(contig):
                        target.write(read)
        else:
            print(f"Streaming unindexed {stem} ({int(bam_bytes) / 1e6:.1f} MB); "
                  "retaining mitochondrial reads only", flush=True)
            stream_unindexed_mt(bam_url, pending, bam_md5, int(bam_bytes))
        pending.replace(subset)
    acquisition_path.write_text(json.dumps(acquisition, indent=2) + "\n", encoding="utf-8")
    if not Path(str(subset) + ".bai").exists():
        pysam.index(str(subset))
    with pysam.AlignmentFile(str(subset), "rb") as source:
        candidates = [name for name in source.references if name in {"MT", "chrM", "M", "chrMT"}]
        if len(candidates) != 1:
            raise ValueError(f"Invalid mitochondrial subset reference: {stem}")
        contig = candidates[0]
        sq = next(row for row in source.header.to_dict()["SQ"] if row["SN"] == contig)
        mapped_reads = sum(1 for _ in source.fetch(contig))
        primary_reads = sum(passing_read(read) for read in source.fetch(contig))
        coverage = source.count_coverage(
            contig, start=0, stop=sq["LN"], quality_threshold=20, read_callback=passing_read,
        )
        depth = [sum(counts) for counts in zip(*coverage, strict=True)]
    return {
        **individual, "study": STUDY, "run_accession": run["run_accession"],
        "sample_accession": run["sample_accession"],
        "secondary_sample_accession": run["secondary_sample_accession"],
        "ena_sample_alias": run["sample_alias"], "status": "downloaded_mt_read_subset",
        "source_bam_url": bam_url, "source_bam_bytes": int(bam_bytes),
        "source_bam_declared_md5": bam_md5, **acquisition,
        "mt_subset_path": str(subset), "mt_subset_bytes": subset.stat().st_size,
        "mt_subset_sha256": digest(subset), "mt_subset_index_sha256": digest(Path(str(subset) + ".bai")),
        "mitochondrial_reference_header": sq,
        "mapped_mt_reads": mapped_reads, "mapq30_primary_nonduplicate_mt_reads": primary_reads,
        "coverage_filters": {
            "minimum_mapping_quality": 30, "minimum_base_quality": 20,
            "excluded_sam_flag_mask": "0xF04", "terminal_damage_trimming": False,
            "overlapping_mates": "read-level depth; not fragment-overlap-corrected",
        },
        "mean_filtered_mt_depth": sum(depth) / len(depth),
        "positions_at_depth": {str(threshold): sum(value >= threshold for value in depth)
                               for threshold in (1, 3, 5, 10, 20)},
        "consensus_generated": False,
        "sequence_status": "Aligned mitochondrial reads only; not a complete consensus mitogenome",
    }


def main(output: Path = DEFAULT_OUTPUT, workers: int = 4) -> None:
    if not 1 <= workers <= 8:
        raise ValueError("Use between one and eight download workers")
    output = output.resolve()
    metadata = output / "source_metadata"
    for directory in (metadata, output / "source_indexes", output / "mt_reads"):
        directory.mkdir(parents=True, exist_ok=True)
    table = metadata / "Wang2021_OnlineTable1.xlsx"
    fetch_table(table)
    ena = metadata / f"ena_{STUDY}.tsv"
    if not ena.exists():
        ena.write_bytes(fetch(ENA_URL))
    individuals = individual_metadata(table)
    (metadata / "Wang2021_Hanben_Gongguan_individuals.json").write_text(
        json.dumps(individuals, indent=2) + "\n", encoding="utf-8",
    )
    runs = list(csv.DictReader(io.StringIO(ena.read_text()), delimiter="\t"))
    tasks: list[tuple[dict[str, Any], dict[str, str]]] = []
    for individual in individuals:
        matching = [run for run in runs if run["sample_alias"] == individual["individual_id"]]
        individual["ena_run_accessions"] = ";".join(run["run_accession"] for run in matching)
        individual["ena_sample_accessions"] = ";".join(run["sample_accession"] for run in matching)
        if not matching:
            raise ValueError(f"No exact ENA individual alias match: {individual['individual_id']}")
        tasks.extend((individual, run) for run in matching)
    with (output / "individuals.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(individuals[0]))
        writer.writeheader()
        writer.writerows(individuals)
    manifest_path = output / "manifest.json"
    results: list[dict[str, Any]] = []
    manifest: dict[str, Any] = {
        "study": STUDY, "paper_url": PAPER, "publication_doi": "10.1038/s41586-021-03336-2",
        "metadata_url": TABLE_URL, "metadata_sha256": digest(table),
        "ena_report_url": ENA_URL, "ena_report_sha256": digest(ena),
        "individual_count": len(individuals), "run_count": len(tasks),
        "site_counts": dict(Counter(row["site"] for row in individuals)),
        "acquisition_method": (
            "Public HTTPS BAI download and indexed mitochondrial range extraction; "
            "unindexed BAMs are stream-filtered without retaining nuclear reads"
        ),
        "consensus_generated": False,
        "complete_mitogenome_fasta_status": "No author-released consensus FASTA identified in this study deposit",
        "limitations": [
            "BAI checksums are verified; full remote nuclear BAM MD5s are declared but not verified.",
            "Reads remain aligned to their deposited reference; identity must be checked before consensus calling.",
            "Reported paper mtDNA depth may exceed depth in deposited mixed nuclear BAMs.",
            "No consensus bases are invented or reference-filled; missing data remain missing.",
            "I13722 appears in secondary databases but is absent from this original primary table and ENA study report.",
        ],
        "runs": results,
    }
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(extract_run, person, run, output): (person, run) for person, run in tasks}
        for future in as_completed(futures):
            person, run = futures[future]
            try:
                result = future.result()
                print(f"Downloaded {person['individual_id']} {run['run_accession']}: "
                      f"{result['mapped_mt_reads']} mt reads", flush=True)
            except Exception as error:
                # Keep a per-run failure ledger while allowing independent public
                # downloads to finish; the process exits unsuccessfully below.
                result = {**person, "run_accession": run["run_accession"],
                          "status": "failed", "error": f"{type(error).__name__}: {error}"}
                print(f"FAILED {person['individual_id']}: {error}", flush=True)
            results.append(result)
            results.sort(key=lambda row: (row["individual_id"], row["run_accession"]))
            manifest["retrieved_at_utc"] = datetime.now(UTC).isoformat()
            manifest["downloaded_run_count"] = sum(row["status"] != "failed" for row in results)
            manifest["failed_run_count"] = sum(row["status"] == "failed" for row in results)
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if manifest["failed_run_count"]:
        raise typer.Exit(1)
    print(f"Complete: {len(results)} mitochondrial read subsets; no consensus inferred.")


if __name__ == "__main__":
    typer.run(main)
