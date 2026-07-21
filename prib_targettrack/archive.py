"""Archive inventory and schema discovery without a guessed TargetTrack XPath."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
import tarfile
import xml.etree.ElementTree as element_tree
from collections import Counter, defaultdict
from pathlib import Path
from typing import BinaryIO


def file_md5(path: str | Path) -> str:
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def archive_inventory(archive_path: str | Path, output_path: str | Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive.getmembers():
            rows.append(
                {
                    "member_name": member.name,
                    "size_bytes": member.size,
                    "is_file": member.isfile(),
                    "is_directory": member.isdir(),
                }
            )
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["member_name"])
        writer.writeheader()
        writer.writerows(rows)
    return rows


def extract_documentation(archive_path: str | Path, output_dir: str | Path) -> list[str]:
    """Extract source documentation only; never expand the raw XML or contacts."""
    wanted_names = {
        "README_tt.txt",
        "targetTrack-v1.4.1.xsd",
        "targetTrackEnumeratedDataItems-v1.4.1.xls",
        "targetTrack-v1.4.1.pdf",
    }
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    extracted_paths: list[str] = []
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile() or Path(member.name).name not in wanted_names:
                continue
            source = archive.extractfile(member)
            if source is None:
                continue
            target = destination / Path(member.name).name
            with source, target.open("wb") as handle:
                shutil.copyfileobj(source, handle)
            extracted_paths.append(str(target))
    return sorted(extracted_paths)


def _open_xml_member(archive: tarfile.TarFile) -> tuple[str, BinaryIO]:
    candidates = [member for member in archive.getmembers() if member.isfile() and member.name.casefold().endswith(".xml.gz")]
    if not candidates:
        raise ValueError("No compressed XML member found in archive inventory.")
    member = next((candidate for candidate in candidates if candidate.name.casefold().endswith("/tt.xml.gz")), candidates[0])
    extracted = archive.extractfile(member)
    if extracted is None:
        raise ValueError(f"Could not read XML member: {member.name}")
    return member.name, extracted


def schema_snapshot(
    archive_path: str | Path,
    output_path: str | Path,
    sample_records: int = 200,
) -> dict[str, object]:
    """Collect names and limited status-like values, never generic free text."""
    tag_counts: Counter[str] = Counter()
    status_values: dict[str, set[str]] = defaultdict(set)
    parsed_elements = 0
    element_limit = max(500, sample_records * 200)
    with tarfile.open(archive_path, "r:gz") as archive:
        member_name, compressed = _open_xml_member(archive)
        with compressed, gzip.GzipFile(fileobj=compressed, mode="rb") as xml_handle:
            for _, element in element_tree.iterparse(xml_handle, events=("end",)):
                name = local_name(element.tag)
                tag_counts[name] += 1
                name_key = name.casefold()
                text = (element.text or "").strip()
                if text and len(text) <= 120 and any(token in name_key for token in ("status", "stage", "stop")):
                    status_values[name].add(text)
                element.clear()
                parsed_elements += 1
                if parsed_elements >= element_limit:
                    break
    snapshot = {
        "archive_member": member_name,
        "sample_records_requested": sample_records,
        "sampled_xml_end_elements": parsed_elements,
        "tag_counts": dict(sorted(tag_counts.items())),
        "status_or_stage_examples": {key: sorted(values)[:25] for key, values in sorted(status_values.items())},
        "privacy_note": "Only short values from status/stage/stop-like fields are retained in this discovery snapshot.",
    }
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(snapshot, indent=2, sort_keys=True), encoding="utf-8")
    return snapshot
