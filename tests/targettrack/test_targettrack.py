from __future__ import annotations

import csv
import gzip
import io
import json
import tarfile
import tempfile
import unittest
from pathlib import Path

from prib_targettrack.analysis import run_analysis, run_streaming_analysis
from prib_targettrack.archive import archive_inventory, extract_documentation, schema_snapshot
from prib_targettrack.statistics import jeffreys_interval, wilson_interval
from prib_targettrack.xml_stream import stream_targettrack_events


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = REPOSITORY_ROOT / "data" / "fixtures" / "targettrack" / "synthetic_events.csv"


class TargetTrackFunnelTests(unittest.TestCase):
    def test_monotonic_funnel_and_stop_boundaries(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = run_analysis(FIXTURE, temporary)
            self.assertEqual(result["trial_units"], 4)
            with (Path(temporary) / "funnel_trial_level.csv").open(encoding="utf-8") as handle:
                funnel = {row["canonical_stage"]: row for row in csv.DictReader(handle)}
            self.assertEqual(funnel["purified"]["observed_passed"], "2")
            self.assertEqual(funnel["expressed"]["monotonic_passed"], "2")
            self.assertEqual(funnel["expressed"]["inferred_prerequisite_passes"], "1")
            with (Path(temporary) / "analysis_units_trial_level.csv").open(encoding="utf-8") as handle:
                units = {row["target_id"]: row for row in csv.DictReader(handle)}
            self.assertEqual(units["T-002"]["technical_failure_code"], "expression_failed")
            self.assertEqual(units["T-004"]["terminal_stop_class"], "unknown_stop")

    def test_retry_and_later_success_are_not_erased_or_counted_as_terminal_failure(self) -> None:
        rows = [
            ["T-005", "A", "C-005A", "2014-05-01", "selected", ""],
            ["T-005", "A", "C-005A", "2014-05-02", "cloned", ""],
            ["T-005", "A", "C-005A", "2014-05-03", "", "expression failed"],
            ["T-005", "B", "C-005B", "2014-05-04", "selected", ""],
            ["T-005", "B", "C-005B", "2014-05-05", "cloned", ""],
            ["T-005", "B", "C-005B", "2014-05-06", "expressed", ""],
            ["T-005", "B", "C-005B", "2014-05-07", "soluble", ""],
            ["T-005", "B", "C-005B", "2014-05-08", "purified", ""],
            ["T-006", "A", "C-006", "2014-06-01", "selected", ""],
            ["T-006", "A", "C-006", "2014-06-02", "cloned", ""],
            ["T-006", "A", "C-006", "2014-06-03", "", "work stopped"],
            ["T-006", "A", "C-006", "2014-06-04", "deposited", ""],
        ]
        with tempfile.TemporaryDirectory() as temporary:
            events_path = Path(temporary) / "events.csv"
            with events_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["target_id", "trial_id", "construct_id", "event_timestamp", "raw_status", "raw_stop_status"])
                writer.writerows(rows)
            run_analysis(events_path, temporary)
            with (Path(temporary) / "analysis_units_target_level.csv").open(encoding="utf-8") as handle:
                targets = {row["target_id"]: row for row in csv.DictReader(handle)}
            self.assertEqual(targets["T-005"]["mixed_outcome_flag"], "True")
            self.assertEqual(targets["T-005"]["highest_monotonic_stage"], "purified")
            with (Path(temporary) / "analysis_units_trial_level.csv").open(encoding="utf-8") as handle:
                units = {row["unit_id"]: row for row in csv.DictReader(handle)}
            self.assertEqual(units["T-006|A|C-006"]["contradiction_flag"], "True")
            self.assertEqual(units["T-006|A|C-006"]["terminal_stop_class"], "not_stopped")

    def test_unknown_status_blocks_strict_analysis_and_is_auditable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            events_path = Path(temporary) / "events.csv"
            events_path.write_text(
                "target_id,trial_id,construct_id,raw_status\nT-007,A,C-007,not-a-reviewed-status\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "Unmapped raw statuses"):
                run_analysis(events_path, temporary)
            self.assertIn("not-a-reviewed-status", (Path(temporary) / "unmapped_statuses.csv").read_text(encoding="utf-8"))

    def test_streaming_aggregation_matches_fixture_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            summary = run_streaming_analysis(FIXTURE, temporary)
            self.assertEqual(summary["execution_mode"], "streaming")
            self.assertEqual(summary["trial_units"], 4)
            self.assertEqual(summary["target_units"], 4)
            with (Path(temporary) / "transition_outcomes.csv").open(encoding="utf-8") as handle:
                transitions = list(csv.DictReader(handle))
            self.assertEqual(transitions[1]["unknown_stops"], "1")
            self.assertTrue((Path(temporary) / "data_quality_report.json").exists())

    def test_intervals_are_bounded(self) -> None:
        low, high = wilson_interval(3, 10)
        mean, beta_low, beta_high = jeffreys_interval(3, 7)
        self.assertTrue(0.0 <= low <= high <= 1.0)
        self.assertTrue(0.0 <= beta_low <= mean <= beta_high <= 1.0)

    def test_archive_inspection_uses_inventory_then_safe_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            archive_path = Path(temporary) / "mini.tar.gz"
            xml_bytes = b"<root><target><status>selected</status><stopStatus>work stopped</stopStatus></target></root>"
            compressed = gzip.compress(xml_bytes)
            with tarfile.open(archive_path, "w:gz") as archive:
                metadata = tarfile.TarInfo("TargetTrack XML files/tt.xml.gz")
                metadata.size = len(compressed)
                archive.addfile(metadata, io.BytesIO(compressed))
            inventory_path = Path(temporary) / "inventory.csv"
            snapshot_path = Path(temporary) / "snapshot.json"
            archive_inventory(archive_path, inventory_path)
            self.assertEqual(extract_documentation(archive_path, Path(temporary) / "docs"), [])
            snapshot = schema_snapshot(archive_path, snapshot_path, sample_records=1)
            self.assertTrue(inventory_path.exists())
            self.assertIn("status", snapshot["status_or_stage_examples"])
            self.assertEqual(json.loads(snapshot_path.read_text(encoding="utf-8"))["archive_member"], "TargetTrack XML files/tt.xml.gz")

    def test_stream_parser_keeps_experimental_history_and_omits_contact_data(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            archive_path = Path(temporary) / "mini.tar.gz"
            xml_bytes = (
                b"<targetTrack><target id='ABC-T1'><targetId>T1</targetId><trialList>"
                b"<trial id='1'><statusHistoryList>"
                b"<statusHistory id='1'><status>selected</status><dateComplete>2017-01-01</dateComplete></statusHistory>"
                b"<statusHistory id='2'><status>work stopped</status><dateComplete>2017-01-02</dateComplete></statusHistory>"
                b"</statusHistoryList><trialSequenceList><trialSequence id='2'><oneLetterCode>MAAA</oneLetterCode></trialSequence>"
                b"</trialSequenceList></trial></trialList></target></targetTrack>"
            )
            compressed = gzip.compress(xml_bytes)
            with tarfile.open(archive_path, "w:gz") as archive:
                metadata = tarfile.TarInfo("TargetTrack XML files/tt.xml.gz")
                metadata.size = len(compressed)
                archive.addfile(metadata, io.BytesIO(compressed))
            events_path = Path(temporary) / "events.csv"
            summary = stream_targettrack_events(archive_path, events_path)
            self.assertEqual(summary["events"], 2)
            text = events_path.read_text(encoding="utf-8")
            self.assertIn("work stopped", text)
            self.assertNotIn("oneLetterCode", text)


if __name__ == "__main__":
    unittest.main()
