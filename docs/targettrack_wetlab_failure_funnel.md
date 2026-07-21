# TargetTrack wet-lab failure funnel

## Purpose

PRIB-KI now includes a CPU-only, reproducible foundation for analysing
historical protein-production and structure-determination progression in the
Protein Structure Initiative's TargetTrack archive. The archive is a starting
evidence layer for descriptive historical priors; it is not treated as
universal industrial wet-lab ground truth.

The source is the final 30 June 2017 TargetTrack snapshot available from
Zenodo: https://doi.org/10.5281/zenodo.821654. Its record describes the XML
archive, documentation, enumeration spreadsheet, published MD5 and CC BY-SA
4.0 licence. Attribution and data-handling rules are recorded in
[THIRD_PARTY_DATA.md](../THIRD_PARTY_DATA.md).

## What is implemented now

- verified downloader with MD5 and immutable provenance output;
- tar inventory and conservative XML schema discovery, before any fixed XPath
  parser is assumed;
- reviewed-map workflow for raw statuses and stop reasons;
- a normalized-event analysis core with separate trial/construct and target
  aggregation;
- observed-only and monotonic funnels, so inferred prerequisites are never
  presented as directly observed data;
- separate technical failure, nontechnical stop, unknown-stop and censored
  accounting;
- adjacent progression rates with Wilson intervals and Jeffreys-prior
  beta-binomial intervals;
- synthetic fixtures and CPU-only tests.

## Local 2017 snapshot run

On 2026-07-21, the official archive was downloaded locally, checked against
its published MD5, inspected with its supplied documentation, and parsed
without retaining contact information or protocol free text. The local,
Git-ignored run contains 3,783,070 status-history events for 961,548
trial/construct units and 335,771 targets. The complete raw-status inventory
was mapped with zero unmapped status or stop-reason values.

This source snapshot does **not** contain an explicit technical or
nontechnical terminal-failure classification in the parsed status history.
It contains 173,989 terminal trial records labelled \`work stopped\`, which
remain \`unknown_stop\` rather than molecular failures. Consequently,
\`resolved_probability\` and beta-posterior fields are deliberately blank in
the local transition output: they require at least one explicitly technical
terminal failure at the transition. The remaining \`raw_probability\` field is
only a descriptive historical progression fraction; its denominator includes
unknown stops and censored records and it is not a calibrated success
probability.

The XML schema, documentation and raw enumeration have now been inspected.
The software map covers the observed values, but its scientific interpretation
still requires biological and statistical review before any training or
comparative performance claim.

## Canonical interpretation

The shared progression backbone is:

    selected -> cloned -> expressed -> soluble -> purified
      -> crystallized -> diffraction_quality -> structure_determined -> deposited

X-ray, NMR and EM branch details remain separate where appropriate. A later
observed stage can imply that prerequisite work happened, but the output marks
such stages as inferred prerequisites. A record saying “work stopped” remains
unknown unless a deterministic, reviewed rule classifies a reason; it is never
silently counted as molecular failure.

The first release will name any output a “TargetTrack historical progression
prior” or “structural-genomics experimental progression probability.” It will
not describe it as an industrial wet-lab success probability.

## Safe local workflow

From the repository root:

    python -m prib_targettrack inspect --archive data/raw/targettrack/TargetTrack-1Jul2017.tar.gz
    python -m prib_targettrack parse --archive data/raw/targettrack/TargetTrack-1Jul2017.tar.gz --out data/interim/targettrack/events.csv
    python -m prib_targettrack analyze --events data/interim/targettrack/events.csv --out reports/targettrack
    python -m prib_targettrack analyze --events data/fixtures/targettrack/synthetic_events.csv

Download is intentionally explicit:

    python -m prib_targettrack download

It downloads to an ignored directory and refuses a bad MD5. The raw archive,
interim records, processed row-level data and TargetTrack-derived reports are
ignored by Git. Public release of substantial transformed row-level data needs
a separate CC BY-SA 4.0 data-release decision and legal review.

## Evidence boundary

TargetTrack records span historical structural-genomics centres, protocols and
years. They are useful for failure taxonomy, censoring-aware funnel design and
historical progression priors. They cannot on their own establish
partner-specific performance, modern manufacturing success, calibrated
candidate probabilities or clinical outcomes. The PRIB-KI demonstrator remains
unchanged: it produces relative screening indices, not TargetTrack-derived
scores.
