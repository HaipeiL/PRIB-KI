"""Configuration loading for the TargetTrack workflow.

The checked-in .yaml file is JSON-formatted YAML, so the scaffold stays
dependency-free before the data workflow is operational.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = REPOSITORY_ROOT / "configs" / "targettrack.yaml"


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    config_path = Path(path) if path else DEFAULT_CONFIG_PATH
    with config_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def repository_path(config: dict[str, Any], key: str) -> Path:
    return REPOSITORY_ROOT / config["paths"][key]
