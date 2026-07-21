"""Command line contract for the CPU-only TargetTrack foundation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .analysis import normalize_events, read_event_rows, run_analysis, write_rows
from .archive import archive_inventory, extract_documentation, schema_snapshot
from .config import load_config, repository_path
from .download import download_targettrack
from .mapping import StageMapper, StopMapper
from .xml_stream import stream_targettrack_events


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(prog="python -m prib_targettrack")
    command.add_argument("--config", help="Path to JSON-compatible YAML configuration.")
    subcommands = command.add_subparsers(dest="command", required=True)
    download = subcommands.add_parser("download", help="Download TargetTrack with MD5 verification.")
    download.add_argument("--out", help="Ignored archive destination.")
    inspect = subcommands.add_parser("inspect", help="Inventory archive and create a schema snapshot.")
    inspect.add_argument("--archive", required=True)
    inspect.add_argument("--sample-records", type=int, default=200)
    inspect.add_argument("--out", help="Report directory.")
    parse = subcommands.add_parser("parse", help="Stream TargetTrack XML to non-contact status-history events.")
    parse.add_argument("--archive", required=True)
    parse.add_argument("--out", help="Ignored interim event CSV path.")
    analyze = subcommands.add_parser("analyze", help="Run the normalized-event funnel core.")
    analyze.add_argument("--events", required=True, help="CSV or JSONL event file using the documented normalized fields.")
    analyze.add_argument("--out", help="Report directory.")
    analyze.add_argument("--exploratory", action="store_true", help="Allow unmapped statuses while writing their review file.")
    audit = subcommands.add_parser("audit-mappings", help="Write status and stop-reason review inventories.")
    audit.add_argument("--events", required=True)
    audit.add_argument("--out", help="Report directory.")
    return command


def main(argv: list[str] | None = None) -> int:
    arguments = parser().parse_args(argv)
    config = load_config(arguments.config)
    report_dir = Path(getattr(arguments, "out", None) or repository_path(config, "report_dir"))
    if arguments.command == "download":
        print(json.dumps(download_targettrack(config, arguments.out), indent=2))
        return 0
    if arguments.command == "inspect":
        report_dir.mkdir(parents=True, exist_ok=True)
        inventory = archive_inventory(arguments.archive, report_dir / "archive_inventory.csv")
        documentation = extract_documentation(arguments.archive, report_dir / "source_documentation")
        snapshot = schema_snapshot(arguments.archive, report_dir / "schema_snapshot.json", arguments.sample_records)
        print(json.dumps({"inventory_members": len(inventory), "documentation": documentation, "snapshot": snapshot}, indent=2))
        return 0
    if arguments.command == "parse":
        output_path = Path(arguments.out) if arguments.out else repository_path(config, "interim_dir") / "events.csv"
        print(json.dumps(stream_targettrack_events(arguments.archive, output_path), indent=2))
        return 0
    if arguments.command == "analyze":
        summary = run_analysis(arguments.events, report_dir, strict=not arguments.exploratory)
        print(json.dumps(summary, indent=2))
        return 0
    rows = read_event_rows(arguments.events)
    _, statuses, stops = normalize_events(rows, StageMapper(), StopMapper())
    report_dir.mkdir(parents=True, exist_ok=True)
    write_rows(
        report_dir / "unmapped_statuses.csv",
        [{"raw_status": key, "count": value} for key, value in sorted(statuses.items())],
        ["raw_status", "count"],
    )
    write_rows(
        report_dir / "unmapped_stop_reasons.csv",
        [{"raw_stop_reason": key, "count": value} for key, value in sorted(stops.items())],
        ["raw_stop_reason", "count"],
    )
    print(json.dumps({"unmapped_statuses": dict(statuses), "unmapped_stops": dict(stops)}, indent=2))
    return 0
