"""Streaming extraction of non-contact TargetTrack experimental history."""

from __future__ import annotations

import csv
import gzip
import json
import os
import tarfile
import tempfile
import xml.etree.ElementTree as element_tree
from collections import Counter
from pathlib import Path


EVENT_FIELDS = [
    "target_id",
    "trial_id",
    "construct_id",
    "center_id",
    "event_id",
    "event_timestamp",
    "raw_status",
    "raw_stop_status",
    "raw_stop_detail",
    "source_xml_path",
]


def _text(element: element_tree.Element, path: str) -> str:
    value = element.findtext(path)
    return value.strip() if value else ""


def _xml_member(archive: tarfile.TarFile) -> tarfile.TarInfo:
    for member in archive.getmembers():
        if member.isfile() and member.name.endswith("/tt.xml.gz"):
            return member
    raise ValueError("TargetTrack tt.xml.gz was not found in the archive.")


def stream_targettrack_events(archive_path: str | Path, output_path: str | Path) -> dict[str, object]:
    """Write trial status-history events without contact or protocol free text."""
    archive_file = Path(archive_path)
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=destination.name + ".", suffix=".part", dir=destination.parent)
    status_counts: Counter[str] = Counter()
    stop_counts: Counter[str] = Counter()
    target_count = 0
    trial_count = 0
    event_count = 0
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=EVENT_FIELDS)
            writer.writeheader()
            with tarfile.open(archive_file, "r:gz") as archive:
                member = _xml_member(archive)
                extracted = archive.extractfile(member)
                if extracted is None:
                    raise ValueError(f"Could not read XML member: {member.name}")
                with extracted, gzip.GzipFile(fileobj=extracted, mode="rb") as xml_handle:
                    root: element_tree.Element | None = None
                    for event, element in element_tree.iterparse(xml_handle, events=("start", "end")):
                        if root is None and event == "start":
                            root = element
                        if event != "end" or element.tag != "target":
                            continue
                        target_count += 1
                        target_id = _text(element, "targetId") or element.attrib.get("id", "")
                        center_id = element.attrib.get("id", "").split("-", 1)[0]
                        for trial_index, trial in enumerate(element.findall("./trialList/trial"), start=1):
                            trial_count += 1
                            trial_id = trial.attrib.get("id") or _text(trial, "trialId") or str(trial_index)
                            sequence_ids = [
                                sequence.attrib.get("id", "")
                                for sequence in trial.findall("./trialSequenceList/trialSequence")
                                if sequence.attrib.get("id")
                            ]
                            construct_id = ";".join(sequence_ids) or trial_id
                            histories = trial.findall("./statusHistoryList/statusHistory")
                            for history_index, history in enumerate(histories, start=1):
                                raw_status = _text(history, "status")
                                raw_stop_status = _text(history, "stopStatus")
                                raw_stop_detail = ""
                                if not raw_status and not raw_stop_status:
                                    continue
                                event_count += 1
                                event_id = history.attrib.get("id") or str(history_index)
                                status_counts[raw_status] += bool(raw_status)
                                stop_counts[raw_stop_status] += bool(raw_stop_status)
                                writer.writerow(
                                    {
                                        "target_id": target_id,
                                        "trial_id": trial_id,
                                        "construct_id": construct_id,
                                        "center_id": center_id,
                                        "event_id": event_id,
                                        "event_timestamp": _text(history, "dateComplete"),
                                        "raw_status": raw_status,
                                        "raw_stop_status": raw_stop_status,
                                        "raw_stop_detail": raw_stop_detail,
                                        "source_xml_path": member.name,
                                    }
                                )
                        element.clear()
                        if root is not None:
                            root.clear()
        os.replace(temporary_name, destination)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)

    inventory_path = destination.with_name("raw_status_inventory.csv")
    with inventory_path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=["field", "raw_value", "count"])
        writer.writeheader()
        for value, count in sorted(status_counts.items()):
            writer.writerow({"field": "status", "raw_value": value, "count": count})
        for value, count in sorted(stop_counts.items()):
            writer.writerow({"field": "stopStatus", "raw_value": value, "count": count})
    summary = {
        "archive": str(archive_file),
        "xml_member": member.name,
        "events": event_count,
        "targets": target_count,
        "trials": trial_count,
        "events_path": str(destination),
        "status_inventory_path": str(inventory_path),
        "privacy_note": "Contact information, protocol text, addresses, names and email fields were not extracted.",
    }
    destination.with_name("parse_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return summary
