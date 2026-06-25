# Roadmap

Maintained by **PRIB-KI Lab**.

This roadmap tracks technical risk reduction rather than feature count.

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
- sequence-clustered train and test splits;
- simple physicochemical baselines;
- benchmark scripts;
- uncertainty and confidence intervals;
- documented failure cases.

Exit criteria:

- at least one experimental endpoint is evaluated on an independent holdout set;
- results are compared with simple baselines;
- data leakage checks are documented.

## v0.4 — calibrated research prototype

Deliverables:

- experimentally calibrated risk modules;
- protein-class-specific applicability checks;
- structure-derived features where justified;
- model and data versioning;
- machine-readable result reports;
- feature ablation and sensitivity analysis.

Exit criteria:

- reproducible improvement over baseline methods;
- stable performance under family-aware splitting;
- calibrated or explicitly non-probabilistic outputs;
- out-of-domain candidates are flagged.

## v0.5 — prospective validation

Deliverables:

- frozen prediction pipeline;
- blinded candidate panel;
- predefined assay package;
- prediction registry before assay execution;
- prospective evaluation report.

Exit criteria:

- predictions are generated before wet-lab results are available;
- all protocol deviations are recorded;
- both successful and failed predictions are analyzed.

## v0.6 — pilot integration

Deliverables:

- batch API;
- project-level configuration;
- audit-ready reports;
- role-based access and data handling rules;
- pilot-specific monitoring;
- assay recommendation output.

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