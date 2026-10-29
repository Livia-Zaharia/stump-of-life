"""Dataset acquisition tests use synthetic local FASTA, never the network."""

import hashlib
import io
import json
from datetime import datetime
from pathlib import Path
from urllib.error import URLError

import pytest

from stump_of_life import datasets
from stump_of_life.io import read_fasta


def panel_payload(*, stop: int = 541056) -> bytes:
    return "".join(
        f">ENA|KF{number}|KF{number}.1 synthetic description\naCgTnR\n"
        for number in range(540505, stop)
    ).encode("ascii")


def test_parse_ena_headers_wrapping_and_iupac_are_preserved() -> None:
    records = datasets.parse_ena_fasta(
        ">ENA|KF540505|KF540505.1 description\nacgtn\nry sw\n"
        ">ENA|KF540506|KF540506.1 another description\nT\n"
    )
    assert records == {"KF540505.1": "ACGTNRYSW", "KF540506.1": "T"}


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ("", "No ENA FASTA"),
        ("ACGT\n", "before an ENA FASTA header"),
        (">KF540505.1 description\nACGT\n", "Invalid ENA FASTA header"),
        (">ENA|KF540505|KF540506.1\nACGT\n", "matching accessions"),
        (">ENA|KF540505|KF540505.1\n", "Empty ENA FASTA sequence"),
        (">ENA|KF540505|KF540505.1\nAC-GT\n", "Invalid ungapped IUPAC"),
        (">ENA|KF540505|KF540505.1\nACXT\n", "Invalid ungapped IUPAC"),
        (
            ">ENA|KF540505|KF540505.1\nA\n>ENA|KF540505|KF540505.1\nC\n",
            "Duplicate ENA accession",
        ),
    ],
)
def test_parse_invalid_payloads_fail_without_repair(payload: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        datasets.parse_ena_fasta(payload)


def test_fetch_complete_panel_saves_exact_response_and_sequence_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload = panel_payload()
    requests: list[tuple[object, int]] = []

    def fake_urlopen(request: object, timeout: int) -> io.BytesIO:
        requests.append((request, timeout))
        return io.BytesIO(payload)

    monkeypatch.setattr(datasets, "urlopen", fake_urlopen)
    output = tmp_path / "panel"
    summary = datasets.fetch_taiwan_panel(output)
    assert len(requests) == 1
    request, timeout = requests[0]
    assert request.full_url == datasets.TAIWAN_PANEL_URL
    assert "stump-of-life" in request.get_header("User-agent")
    assert timeout == 60
    assert (output / "ena_original_551.fasta").read_bytes() == payload
    records = read_fasta(output / "mitogenomes.fasta")
    assert len(records) == 551
    assert list(records)[0] == "KF540505.1"
    assert list(records)[-1] == "KF541055.1"
    assert set(records.values()) == {"ACGTNR"}
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["raw_sha256"] == hashlib.sha256(payload).hexdigest()
    assert provenance["source"] == datasets.TAIWAN_PANEL_URL
    assert provenance["doi"] == "10.1016/j.ajhg.2014.02.003"
    assert provenance["sequences"]["KF540505.1"] == {
        "length": 6,
        "sha256": hashlib.sha256(b"ACGTNR").hexdigest(),
    }
    assert datetime.fromisoformat(provenance["retrieved_utc"]).utcoffset().total_seconds() == 0
    assert provenance["sequence_filling"] is False
    assert summary == {"count": 551, "output": str(output), "sha256": provenance["raw_sha256"]}


@pytest.mark.parametrize(
    "payload",
    [
        panel_payload(stop=541055),
        panel_payload().replace(b"KF541055.1", b"KF541055.2"),
        panel_payload().replace(b"KF541055", b"KF541056"),
    ],
)
def test_fetch_rejects_truncated_or_wrong_version_ranges_before_writing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    payload: bytes,
) -> None:
    monkeypatch.setattr(datasets, "urlopen", lambda request, timeout: io.BytesIO(payload))
    output = tmp_path / "incomplete"
    with pytest.raises(ValueError, match="Expected exactly 551 distinct accessions"):
        datasets.fetch_taiwan_panel(output)
    assert not output.exists()


def test_network_errors_propagate_without_creating_dataset(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail(request: object, timeout: int) -> io.BytesIO:
        raise URLError("synthetic network failure")

    monkeypatch.setattr(datasets, "urlopen", fail)
    with pytest.raises(URLError, match="synthetic network failure"):
        datasets.fetch_taiwan_panel(tmp_path / "failed")
    assert not (tmp_path / "failed").exists()
