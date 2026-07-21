# Roadmap

Maintained by **PRIB-KI Lab**.

This roadmap tracks technical risk reduction rather than feature count. The company-level vision is described in [`long_term_vision.md`](long_term_vision.md).

The AI route is described in detail in [`ai_framework.md`](ai_framework.md). In short, PRIB-KI should evolve from a transparent sequence-only demonstrator into two connected AI systems:

```text
A. Protein physical reliability landscape
   build a map of physical danger zones in AI protein space

B. Lab-in-the-loop physical AI autopilot
   use wet-lab feedback to navigate and update that map
```

## Technical route summary

| Phase | AI objective | Evidence gate |
|---|---|---|
| 0. Transparent demonstrator | Inspectable sequence-only workflow and deterministic rankings | Clean checkout reproduces fixed demo outputs |
| 1. Benchmark foundation | Endpoint-specific datasets, leakage-aware splits, and simple baselines | At least one measured endpoint is tested on an independent holdout set |
| 2. Physical reliability landscape | Add protein language-model embeddings, structure/condition augmentation, and failure anchors | AI landscape improves endpoint signal over physicochemical baselines |
| 3. Danger-zone risk models | Learn physical-risk directions, calibrated thresholds, uncertainty, and applicability-domain flags | Model passes ablation, calibration, and out-of-domain checks |
| 4. Lab-in-the-loop validation | Frozen predictions before wet-lab assays, then calibration update | High-risk ranking enriches predefined experimental failures |
| 5. Pilot integration | API, audit-ready reports, evidence registry, partner configuration, and automation connectors | Partner panel can be processed and traced without manual code changes |

## North-star metrics

| Metric | Why it matters |
|---|---|
| Top-k enrichment of experimental failures | Tests whether PRIB-KI improves candidate prioritization |
| Improvement over physicochemical baselines | Proves the AI landscape adds signal beyond simple descriptors |
| Calibration and uncertainty quality | Prevents overconfident life-science AI claims |
| Applicability-domain detection | Flags candidates outside validated evidence |
| Candidate-to-report turnaround | Measures workflow value for partners |
| Traceability coverage | Ensures every decision links to data, model version, configuration, and assay context |

## v0.2 — reproducible demo

Deliverables:

- stable command and UI workflow;
- fixed random seeds;
- cleaned dependency specification;
- versioned configuration;
- input and output schema documentation;
- basic unit tests;
- demo dataset provenance.

Exit criteria:

- a clean checkout can reproduce the same ranked output;
- all automated tests pass;
- every result records code and configuration versions.

## v0.3 — benchmark build

Deliverables:

- endpoint-specific public datasets;
- benchmark schema for sequence, protein class, endpoint, assay condition, unit, source, and usage rights;
- sequence-clustered train and test splits;
- simple physicochemical baselines;
- benchmark scripts;
- uncertainty and confidence intervals;
- documented failure cases.

Exit criteria:

- at least one experimental endpoint is evaluated on an independent holdout set;
- results are compared with simple baselines;
- data leakage checks are documented.

## v0.4 — physical reliability landscape

Deliverables:

- protein language-model embedding pipeline where justified by benchmarks;
- structure-derived features where structure quality and usage rights allow;
- condition metadata representation for pH, temperature, concentration, buffer, storage, and stress context;
- stable and failed anchor retrieval in embedding space;
- learned risk directions for aggregation, stability, expression, SEC, solubility, or self-interaction endpoints;
- calibrated danger-zone rules;
- protein-class-specific applicability checks;
- uncertainty reporting and out-of-domain flags;
- model and data versioning;
- machine-readable result reports;
- feature ablation and sensitivity analysis.

Exit criteria:

- reproducible improvement over baseline methods;
- stable performance under family-aware splitting;
- calibrated or explicitly non-probabilistic outputs;
- out-of-domain candidates are flagged.

## v0.5 — lab-in-the-loop validation

Deliverables:

- frozen prediction pipeline;
- blinded candidate panel;
- predefined assay package;
- prediction registry before assay execution;
- prospective evaluation report.
- active-learning proposal for the next candidate panel;
- calibration update after assay readout ingestion;
- error taxonomy for missed failures and overcalled risks.

Exit criteria:

- predictions are generated before wet-lab results are available;
- all protocol deviations are recorded;
- both successful and failed predictions are analyzed.
- high-risk groups enrich predefined experimental failures relative to baselines.

## v0.6 — pilot integration

Deliverables:

- batch API;
- project-level configuration;
- audit-ready reports;
- evidence registry for frozen predictions and assay readouts;
- assay recommendation and evidence-gap output;
- ELN/LIMS/export integration;
- robot-ready or CRO-ready assay manifest where partner infrastructure supports it;
- role-based access and data handling rules;
- pilot-specific monitoring;
- partner feedback loop for calibrated module updates.

Exit criteria:

- a partner dataset can be processed without manual code changes;
- outputs can be traced to input, model, and configuration versions;
- pilot users can interpret and review the reported risk drivers.

## Deferred work

The following items are intentionally deferred until benchmark evidence supports them:

- complex deep-learning architectures;
- universal cross-protein-class scoring;
- absolute probability of industrial failure;
- regulatory or clinical decision support;
- automated replacement of wet-lab developability testing.
