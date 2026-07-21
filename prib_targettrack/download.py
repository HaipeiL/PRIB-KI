"""Explicit, checksummed downloader for the ignored TargetTrack source archive."""

from __future__ import annotations

import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

from .archive import file_md5
from .config import repository_path


def download_targettrack(config: dict, destination: str | Path | None = None) -> dict[str, str]:
    dataset = config["dataset"]
    raw_dir = repository_path(config, "raw_dir")
    target = Path(destination) if destination else raw_dir / dataset["archive_name"]
    target.parent.mkdir(parents=True, exist_ok=True)
    expected_md5 = dataset["expected_md5"].casefold()

    if target.exists() and file_md5(target) == expected_md5:
        actual_md5 = expected_md5
    else:
        request = Request(dataset["download_url"], headers={"User-Agent": "PRIB-KI-TargetTrack/0.1"})
        file_descriptor, temporary_name = tempfile.mkstemp(prefix=target.name + ".", suffix=".part", dir=target.parent)
        try:
            with os.fdopen(file_descriptor, "wb") as output, urlopen(request, timeout=90) as response:
                for chunk in iter(lambda: response.read(1024 * 1024), b""):
                    output.write(chunk)
            actual_md5 = file_md5(temporary_name)
            if actual_md5 != expected_md5:
                raise ValueError(f"Checksum mismatch: expected {expected_md5}, received {actual_md5}")
            os.replace(temporary_name, target)
        finally:
            if os.path.exists(temporary_name):
                os.unlink(temporary_name)

    provenance = {
        "dataset": dataset["name"],
        "doi": dataset["doi"],
        "archive_name": target.name,
        "source_url": dataset["download_url"],
        "expected_md5": expected_md5,
        "actual_md5": actual_md5,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "license": dataset["license"],
    }
    provenance_path = target.with_suffix(target.suffix + ".provenance.json")
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8")
    return {"archive": str(target), "provenance": str(provenance_path), "md5": actual_md5}
