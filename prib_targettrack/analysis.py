"""Censoring-aware historical funnel analysis for normalized TargetTrack events."""

from __future__ import annotations

import csv
import json
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from .mapping import StageMapper, StopMapper
from .stages import CANONICAL_STAGES, STAGE_RANK, stage_at, transition_name
from .statistics import jeffreys_interval, wilson_interval


def read_event_rows(path: str | Path) -> list[dict[str, str]]:
    source = Path(path)
    if source.suffix.casefold() == ".jsonl":
        with source.open("r", encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(
    path: str | Path,
    rows: Iterable[dict[str, object]],
    fieldnames: list[str] | None = None,
) -> None:
    rows = list(rows)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = fieldnames or list(dict.fromkeys(field for row in rows for field in row))
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def normalize_events(rows: list[dict[str, str]], stage_mapper: StageMapper, stop_mapper: StopMapper):
    normalized: list[dict[str, object]] = []
    unmapped_statuses: dict[str, int] = defaultdict(int)
    unmapped_stops: dict[str, int] = defaultdict(int)
    for source_order, row in enumerate(rows):
        raw_status = str(row.get("raw_status") or row.get("status") or "").strip()
        raw_stop_status = str(row.get("raw_stop_status") or "").strip()
        raw_stop_detail = str(row.get("raw_stop_detail") or "").strip()
        stop = stop_mapper.lookup(raw_stop_status, raw_stop_detail) if (raw_stop_status or raw_stop_detail) else None
        if stop is None and raw_status:
            stop = stop_mapper.lookup(raw_status, raw_stop_detail)
        stage = stage_mapper.lookup(raw_status) if raw_status else None
        if raw_status and stage is None and stop is None:
            unmapped_statuses[raw_status] += 1
        if (raw_stop_status or raw_stop_detail) and stop is None:
            unmapped_stops[raw_stop_status or raw_stop_detail] += 1
        normalized.append(
            {
                "target_id": str(row.get("target_id") or "").strip(),
                "trial_id": str(row.get("trial_id") or "").strip(),
                "construct_id": str(row.get("construct_id") or "").strip(),
                "event_id": str(row.get("event_id") or "").strip() or str(source_order + 1),
                "event_timestamp": str(row.get("event_timestamp") or "").strip(),
                "source_order": source_order,
                "raw_status": raw_status,
                "canonical_stage": stage.canonical_stage if stage else "",
                "canonical_rank": stage.canonical_rank if stage and stage.canonical_rank is not None else "",
                "branch": stage.branch if stage else "",
                "stage_mapping_confidence": stage.confidence if stage else "",
                "raw_stop_status": raw_stop_status,
                "raw_stop_detail": raw_stop_detail,
                "stop_class": stop.stop_class if stop else ("unknown_stop" if (raw_stop_status or raw_stop_detail) else "not_stopped"),
                "failure_code": stop.failure_code if stop else ("other_unspecified" if (raw_stop_status or raw_stop_detail) else ""),
                "failed_at_stage": stop.failed_at_stage if stop else "",
                "stop_mapping_confidence": stop.confidence if stop else "",
            }
        )
    return normalized, unmapped_statuses, unmapped_stops


def _unit_key(event: dict[str, object]) -> tuple[str, str, str]:
    target = str(event["target_id"] or "").strip()
    if not target:
        raise ValueError("Every analysis event must include target_id.")
    return (target, str(event["trial_id"] or "_"), str(event["construct_id"] or "_"))


def _rank(event: dict[str, object]) -> int | None:
    value = event.get("canonical_rank", "")
    return int(value) if value != "" else None


def _failure_start_rank(event: dict[str, object], observed_before: int) -> int:
    failed_at = str(event.get("failed_at_stage", ""))
    if "_to_" in failed_at:
        start = failed_at.split("_to_", 1)[0]
        if start in STAGE_RANK:
            return STAGE_RANK[start]
    rank = _rank(event)
    return rank if rank is not None else observed_before


def build_trial_units(events: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for event in events:
        grouped[_unit_key(event)].append(event)
    units: list[dict[str, object]] = []
    for (target_id, trial_id, construct_id), history in sorted(grouped.items()):
        history.sort(key=lambda event: (str(event["event_timestamp"]), int(event["source_order"])))
        observed_ranks = {_rank(event) for event in history if _rank(event) is not None}
        highest_rank = max(observed_ranks, default=-1)
        contradiction = False
        terminal_stop: dict[str, object] | None = None
        running_rank = -1
        for index, event in enumerate(history):
            event_rank = _rank(event)
            if event_rank is not None:
                running_rank = max(running_rank, event_rank)
            if event["stop_class"] == "not_stopped":
                continue
            start_rank = _failure_start_rank(event, running_rank)
            later_ranks = [_rank(later) for later in history[index + 1 :] if _rank(later) is not None]
            if later_ranks and max(later_ranks) > start_rank:
                contradiction = True
                continue
            terminal_stop = event
        terminal_class = str(terminal_stop["stop_class"]) if terminal_stop else "not_stopped"
        units.append(
            {
                "unit_id": "|".join((target_id, trial_id, construct_id)),
                "target_id": target_id,
                "trial_id": trial_id,
                "construct_id": construct_id,
                "event_count": len(history),
                "highest_observed_stage": stage_at(highest_rank) if highest_rank >= 0 else "",
                "highest_observed_rank": highest_rank,
                "highest_monotonic_stage": stage_at(highest_rank) if highest_rank >= 0 else "",
                "highest_monotonic_rank": highest_rank,
                "observed_ranks": ";".join(str(rank) for rank in sorted(observed_ranks)),
                "explicit_technical_failure": terminal_class == "technical_failure",
                "technical_failure_code": terminal_stop["failure_code"] if terminal_class == "technical_failure" else "",
                "technical_failure_stage": terminal_stop["failed_at_stage"] if terminal_class == "technical_failure" else "",
                "nontechnical_stop": terminal_class == "nontechnical_stop",
                "unknown_stop": terminal_class == "unknown_stop",
                "work_stopped": terminal_class != "not_stopped",
                "terminal_stop_class": terminal_class,
                "contradiction_flag": contradiction,
                "missing_timestamp_flag": any(not event["event_timestamp"] for event in history),
                "inferred_stage_count": highest_rank + 1 - len(observed_ranks) if highest_rank >= 0 else 0,
            }
        )
    return units


def build_target_units(trial_units: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for unit in trial_units:
        grouped[str(unit["target_id"])].append(unit)
    targets: list[dict[str, object]] = []
    for target_id, units in sorted(grouped.items()):
        representative = max(units, key=lambda unit: int(unit["highest_monotonic_rank"]))
        observed = sorted({int(rank) for unit in units for rank in str(unit["observed_ranks"]).split(";") if rank})
        has_technical_failure = any(bool(unit["explicit_technical_failure"]) for unit in units)
        has_later_success = int(representative["highest_monotonic_rank"]) > 1
        targets.append(
            {
                "unit_id": target_id,
                "target_id": target_id,
                "trial_count": len(units),
                "event_count": sum(int(unit["event_count"]) for unit in units),
                "highest_observed_stage": representative["highest_observed_stage"],
                "highest_observed_rank": representative["highest_observed_rank"],
                "highest_monotonic_stage": representative["highest_monotonic_stage"],
                "highest_monotonic_rank": representative["highest_monotonic_rank"],
                "observed_ranks": ";".join(str(rank) for rank in observed),
                "explicit_technical_failure": has_technical_failure and not has_later_success,
                "technical_failure_code": representative["technical_failure_code"] if not has_later_success else "",
                "technical_failure_stage": representative["technical_failure_stage"] if not has_later_success else "",
                "nontechnical_stop": bool(representative["nontechnical_stop"]) and not has_later_success,
                "unknown_stop": bool(representative["unknown_stop"]) and not has_later_success,
                "work_stopped": bool(representative["work_stopped"]) and not has_later_success,
                "terminal_stop_class": representative["terminal_stop_class"] if not has_later_success else "not_stopped",
                "contradiction_flag": any(bool(unit["contradiction_flag"]) for unit in units),
                "mixed_outcome_flag": has_technical_failure and has_later_success,
                "inferred_stage_count": int(representative["highest_monotonic_rank"]) + 1 - len(observed),
            }
        )
    return targets


def funnel_rows(units: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for rank, stage in enumerate(CANONICAL_STAGES):
        observed = sum(rank in {int(value) for value in str(unit["observed_ranks"]).split(";") if value} for unit in units)
        monotonic = sum(int(unit["highest_monotonic_rank"]) >= rank for unit in units)
        rows.append(
            {
                "canonical_stage": stage,
                "canonical_rank": rank,
                "observed_passed": observed,
                "monotonic_passed": monotonic,
                "inferred_prerequisite_passes": monotonic - observed,
            }
        )
    return rows


def transition_rows(units: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for rank in range(len(CANONICAL_STAGES) - 1):
        transition = transition_name(rank)
        at_risk_units = [unit for unit in units if int(unit["highest_monotonic_rank"]) >= rank]
        advanced = sum(int(unit["highest_monotonic_rank"]) >= rank + 1 for unit in at_risk_units)
        technical = sum(
            unit["terminal_stop_class"] == "technical_failure" and unit["technical_failure_stage"] == transition
            for unit in at_risk_units
        )
        nontechnical = sum(unit["terminal_stop_class"] == "nontechnical_stop" for unit in at_risk_units)
        unknown = sum(unit["terminal_stop_class"] == "unknown_stop" for unit in at_risk_units)
        at_risk = len(at_risk_units)
        censored = max(0, at_risk - advanced - technical - nontechnical - unknown)
        raw_low, raw_high = wilson_interval(advanced, at_risk)
        resolved_total = advanced + technical
        resolved_low, resolved_high = wilson_interval(advanced, resolved_total)
        beta_mean, beta_low, beta_high = jeffreys_interval(advanced, technical)
        rows.append(
            {
                "from_stage": stage_at(rank),
                "to_stage": stage_at(rank + 1),
                "at_risk": at_risk,
                "advanced": advanced,
                "technical_failures": technical,
                "nontechnical_stops": nontechnical,
                "unknown_stops": unknown,
                "censored": censored,
                "unknown_or_censored": unknown + censored,
                "raw_probability": advanced / at_risk if at_risk else "",
                "raw_ci_low": raw_low,
                "raw_ci_high": raw_high,
                "resolved_probability": (
                    advanced / resolved_total if technical and resolved_total else ""
                ),
                "resolved_ci_low": resolved_low if technical else "",
                "resolved_ci_high": resolved_high if technical else "",
                "beta_posterior_mean": beta_mean if technical else "",
                "beta_ci_low": beta_low if technical else "",
                "beta_ci_high": beta_high if technical else "",
            }
        )
    return rows


class _SummaryAccumulator:
    def __init__(self) -> None:
        stage_count = len(CANONICAL_STAGES)
        transition_count = stage_count - 1
        self.observed = [0] * stage_count
        self.monotonic = [0] * stage_count
        self.at_risk = [0] * transition_count
        self.advanced = [0] * transition_count
        self.technical = [0] * transition_count
        self.nontechnical = [0] * transition_count
        self.unknown = [0] * transition_count
        self.censored = [0] * transition_count
        self.terminal_stop_counts: dict[str, int] = defaultdict(int)

    def add(self, unit: dict[str, object]) -> None:
        highest_rank = int(unit["highest_monotonic_rank"])
        observed_ranks = {int(value) for value in str(unit["observed_ranks"]).split(";") if value}
        self.terminal_stop_counts[str(unit["terminal_stop_class"])] += 1
        for rank in range(len(CANONICAL_STAGES)):
            if rank in observed_ranks:
                self.observed[rank] += 1
            if highest_rank >= rank:
                self.monotonic[rank] += 1
        for rank in range(len(CANONICAL_STAGES) - 1):
            if highest_rank < rank:
                continue
            self.at_risk[rank] += 1
            if highest_rank >= rank + 1:
                self.advanced[rank] += 1
            elif unit["terminal_stop_class"] == "technical_failure" and unit["technical_failure_stage"] == transition_name(rank):
                self.technical[rank] += 1
            elif unit["terminal_stop_class"] == "nontechnical_stop":
                self.nontechnical[rank] += 1
            elif unit["terminal_stop_class"] == "unknown_stop":
                self.unknown[rank] += 1
            else:
                self.censored[rank] += 1

    def funnel(self) -> list[dict[str, object]]:
        return [
            {
                "canonical_stage": stage,
                "canonical_rank": rank,
                "observed_passed": self.observed[rank],
                "monotonic_passed": self.monotonic[rank],
                "inferred_prerequisite_passes": self.monotonic[rank] - self.observed[rank],
            }
            for rank, stage in enumerate(CANONICAL_STAGES)
        ]

    def transitions(self) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for rank in range(len(CANONICAL_STAGES) - 1):
            at_risk = self.at_risk[rank]
            advanced = self.advanced[rank]
            technical = self.technical[rank]
            resolved_total = advanced + technical
            raw_low, raw_high = wilson_interval(advanced, at_risk)
            resolved_low, resolved_high = wilson_interval(advanced, resolved_total)
            beta_mean, beta_low, beta_high = jeffreys_interval(advanced, technical)
            rows.append(
                {
                    "from_stage": stage_at(rank),
                    "to_stage": stage_at(rank + 1),
                    "at_risk": at_risk,
                    "advanced": advanced,
                    "technical_failures": technical,
                    "nontechnical_stops": self.nontechnical[rank],
                    "unknown_stops": self.unknown[rank],
                    "censored": self.censored[rank],
                    "unknown_or_censored": self.unknown[rank] + self.censored[rank],
                    "raw_probability": advanced / at_risk if at_risk else "",
                    "raw_ci_low": raw_low,
                    "raw_ci_high": raw_high,
                    "resolved_probability": (
                        advanced / resolved_total if technical and resolved_total else ""
                    ),
                    "resolved_ci_low": resolved_low if technical else "",
                    "resolved_ci_high": resolved_high if technical else "",
                    "beta_posterior_mean": beta_mean if technical else "",
                    "beta_ci_low": beta_low if technical else "",
                    "beta_ci_high": beta_high if technical else "",
                }
            )
        return rows


def _open_atomic_csv(path: Path, fieldnames: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".part", dir=path.parent)
    handle = os.fdopen(descriptor, "w", encoding="utf-8", newline="")
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    return handle, writer, temporary_name


def _write_atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".part", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary_name, path)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)


def run_streaming_analysis(events_path: str | Path, output_dir: str | Path, strict: bool = True) -> dict[str, object]:
    """Aggregate XML-order event CSVs with bounded memory.

    The TargetTrack stream extractor orders rows by target, trial and status
    history. This function validates that ordering through a one-unit buffer,
    writes units atomically, and retains only the current target's trials.
    """
    source = Path(events_path)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    stage_mapper, stop_mapper = StageMapper(), StopMapper()
    unmapped_statuses: dict[str, int] = defaultdict(int)
    unmapped_stops: dict[str, int] = defaultdict(int)
    trial_accumulator, target_accumulator = _SummaryAccumulator(), _SummaryAccumulator()
    event_count = 0
    trial_count = 0
    target_count = 0
    trial_path = output / "analysis_units_trial_level.csv"
    target_path = output / "analysis_units_target_level.csv"
    trial_handle, trial_writer, trial_temporary = _open_atomic_csv(
        trial_path,
        [
            "unit_id", "target_id", "trial_id", "construct_id", "event_count",
            "highest_observed_stage", "highest_observed_rank", "highest_monotonic_stage",
            "highest_monotonic_rank", "observed_ranks", "explicit_technical_failure",
            "technical_failure_code", "technical_failure_stage", "nontechnical_stop",
            "unknown_stop", "work_stopped", "terminal_stop_class", "contradiction_flag",
            "missing_timestamp_flag", "inferred_stage_count",
        ],
    )
    target_handle, target_writer, target_temporary = _open_atomic_csv(
        target_path,
        [
            "unit_id", "target_id", "trial_count", "event_count", "highest_observed_stage",
            "highest_observed_rank", "highest_monotonic_stage", "highest_monotonic_rank",
            "observed_ranks", "explicit_technical_failure", "technical_failure_code",
            "technical_failure_stage", "nontechnical_stop", "unknown_stop", "work_stopped",
            "terminal_stop_class", "contradiction_flag", "mixed_outcome_flag", "inferred_stage_count",
        ],
    )
    current_key: tuple[str, str, str] | None = None
    current_history: list[dict[str, object]] = []
    current_target = ""
    current_target_trials: list[dict[str, object]] = []

    def flush_target() -> None:
        nonlocal target_count, current_target_trials
        if not current_target_trials:
            return
        target_unit = build_target_units(current_target_trials)[0]
        target_writer.writerow(target_unit)
        target_accumulator.add(target_unit)
        target_count += 1
        current_target_trials = []

    def flush_trial() -> None:
        nonlocal trial_count, current_history, current_target
        if not current_history:
            return
        trial_unit = build_trial_units(current_history)[0]
        if current_target and trial_unit["target_id"] != current_target:
            flush_target()
        current_target = str(trial_unit["target_id"])
        trial_writer.writerow(trial_unit)
        trial_accumulator.add(trial_unit)
        current_target_trials.append(trial_unit)
        trial_count += 1
        current_history = []

    try:
        with source.open("r", encoding="utf-8-sig", newline="") as handle:
            for source_order, row in enumerate(csv.DictReader(handle)):
                normalized, statuses, stops = normalize_events([row], stage_mapper, stop_mapper)
                event = normalized[0]
                event["source_order"] = source_order
                for value, count in statuses.items():
                    unmapped_statuses[value] += count
                for value, count in stops.items():
                    unmapped_stops[value] += count
                key = _unit_key(event)
                if current_key is not None and key != current_key:
                    flush_trial()
                current_key = key
                current_history.append(event)
                event_count += 1
        flush_trial()
        flush_target()
        trial_handle.close()
        target_handle.close()
        os.replace(trial_temporary, trial_path)
        os.replace(target_temporary, target_path)
    finally:
        if not trial_handle.closed:
            trial_handle.close()
        if not target_handle.closed:
            target_handle.close()
        if os.path.exists(trial_temporary):
            os.unlink(trial_temporary)
        if os.path.exists(target_temporary):
            os.unlink(target_temporary)

    write_rows(
        output / "unmapped_statuses.csv",
        [{"raw_status": value, "count": count} for value, count in sorted(unmapped_statuses.items())],
        ["raw_status", "count"],
    )
    write_rows(
        output / "unmapped_stop_reasons.csv",
        [{"raw_stop_reason": value, "count": count} for value, count in sorted(unmapped_stops.items())],
        ["raw_stop_reason", "count"],
    )
    if strict and unmapped_statuses:
        raise ValueError("Unmapped raw statuses found. Review configs/targettrack_stage_map.csv before analysis release.")
    trial_funnel, target_funnel = trial_accumulator.funnel(), target_accumulator.funnel()
    transitions = trial_accumulator.transitions()
    write_rows(output / "funnel_trial_level.csv", trial_funnel)
    write_rows(output / "funnel_target_level.csv", target_funnel)
    write_rows(output / "transition_outcomes.csv", transitions)
    write_rows(output / "progression_probabilities_global.csv", transitions)
    _write_atomic_json(
        output / "data_quality_report.json",
        {
            "events": event_count,
            "trial_units": trial_count,
            "target_units": target_count,
            "trial_terminal_stop_classes": dict(sorted(trial_accumulator.terminal_stop_counts.items())),
            "target_terminal_stop_classes": dict(sorted(target_accumulator.terminal_stop_counts.items())),
            "unmapped_statuses": dict(sorted(unmapped_statuses.items())),
            "unmapped_stop_statuses": dict(sorted(unmapped_stops.items())),
            "resolved_probability_policy": (
                "Reported only where at least one explicitly technical terminal stop "
                "is available for the transition; otherwise blank."
            ),
        },
    )
    return {
        "events": event_count,
        "trial_units": trial_count,
        "target_units": target_count,
        "trial_terminal_stop_classes": dict(sorted(trial_accumulator.terminal_stop_counts.items())),
        "unmapped_statuses": dict(unmapped_statuses),
        "unmapped_stops": dict(unmapped_stops),
        "output_dir": str(output),
        "execution_mode": "streaming",
    }


def run_analysis(events_path: str | Path, output_dir: str | Path, strict: bool = True) -> dict[str, object]:
    if Path(events_path).stat().st_size > 64 * 1024 * 1024:
        return run_streaming_analysis(events_path, output_dir, strict)
    raw_rows = read_event_rows(events_path)
    normalized, unmapped_statuses, unmapped_stops = normalize_events(raw_rows, StageMapper(), StopMapper())
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    write_rows(output / "normalized_events.csv", normalized)
    write_rows(
        output / "unmapped_statuses.csv",
        [{"raw_status": key, "count": value} for key, value in sorted(unmapped_statuses.items())],
        ["raw_status", "count"],
    )
    write_rows(
        output / "unmapped_stop_reasons.csv",
        [{"raw_stop_reason": key, "count": value} for key, value in sorted(unmapped_stops.items())],
        ["raw_stop_reason", "count"],
    )
    if strict and unmapped_statuses:
        raise ValueError("Unmapped raw statuses found. Review configs/targettrack_stage_map.csv before analysis release.")
    trial_units = build_trial_units(normalized)
    target_units = build_target_units(trial_units)
    trial_funnel, target_funnel = funnel_rows(trial_units), funnel_rows(target_units)
    transitions = transition_rows(trial_units)
    write_rows(output / "analysis_units_trial_level.csv", trial_units)
    write_rows(output / "analysis_units_target_level.csv", target_units)
    write_rows(output / "funnel_trial_level.csv", trial_funnel)
    write_rows(output / "funnel_target_level.csv", target_funnel)
    write_rows(output / "transition_outcomes.csv", transitions)
    write_rows(output / "progression_probabilities_global.csv", transitions)
    return {
        "events": len(normalized),
        "trial_units": len(trial_units),
        "target_units": len(target_units),
        "unmapped_statuses": dict(unmapped_statuses),
        "unmapped_stops": dict(unmapped_stops),
        "output_dir": str(output),
    }
