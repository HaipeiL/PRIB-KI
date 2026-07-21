"""Deterministic mappings for canonical stages and stop classes."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from .config import REPOSITORY_ROOT
from .stages import STAGE_RANK


def normalize_value(value: object) -> str:
    text = str(value or "").strip().casefold()
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", text)).strip("_")


@dataclass(frozen=True)
class StageMapping:
    canonical_stage: str
    canonical_rank: int | None
    branch: str
    confidence: str


@dataclass(frozen=True)
class StopMapping:
    stop_class: str
    failure_code: str
    failed_at_stage: str
    confidence: str


class StageMapper:
    def __init__(self, path: str | Path | None = None) -> None:
        map_path = Path(path) if path else REPOSITORY_ROOT / "configs" / "targettrack_stage_map.csv"
        self._mappings: dict[str, StageMapping] = {}
        with map_path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                key = row.get("normalized_status") or normalize_value(row.get("raw_status"))
                stage = row.get("canonical_stage", "")
                rank = row.get("canonical_rank", "")
                if stage and stage not in STAGE_RANK:
                    raise ValueError(f"Unknown canonical stage in mapping: {stage}")
                self._mappings[key] = StageMapping(
                    canonical_stage=stage,
                    canonical_rank=int(rank) if rank else None,
                    branch=row.get("branch", "common"),
                    confidence=row.get("mapping_confidence", "unknown"),
                )

    def lookup(self, raw_status: object) -> StageMapping | None:
        return self._mappings.get(normalize_value(raw_status))


class StopMapper:
    def __init__(self, path: str | Path | None = None) -> None:
        map_path = Path(path) if path else REPOSITORY_ROOT / "configs" / "targettrack_stop_map.csv"
        self._exact: dict[str, StopMapping] = {}
        self._patterns: list[tuple[re.Pattern[str], StopMapping]] = []
        with map_path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                mapping = StopMapping(
                    stop_class=row["stop_class"],
                    failure_code=row.get("failure_code", ""),
                    failed_at_stage=row.get("failed_at_stage", ""),
                    confidence=row.get("classification_confidence", "unknown"),
                )
                raw_status = normalize_value(row.get("raw_stop_status"))
                if raw_status:
                    self._exact[raw_status] = mapping
                pattern = (row.get("raw_stop_detail_pattern") or "").strip()
                if pattern:
                    self._patterns.append((re.compile(pattern, re.IGNORECASE), mapping))

    def lookup(self, raw_status: object, detail: object = "") -> StopMapping | None:
        exact = self._exact.get(normalize_value(raw_status))
        if exact:
            return exact
        text = str(detail or "")
        for pattern, mapping in self._patterns:
            if pattern.search(text):
                return mapping
        return None
