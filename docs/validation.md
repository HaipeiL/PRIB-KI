# Validation

Maintained by **PRIB-KI Lab**.

This document defines how PRIB-KI distinguishes completed research evidence from stronger generalization, prospective, and product-validation claims.

## 1. Current evidence status — September 2026

PRIB-KI now has two different validation layers that should not be conflated.

### Public demonstrator

The public Streamlit scoring layer remains a transparent workflow demonstrator based on sequence-derived descriptors, motif/liability proxies, and AI-weighted pseudo-target mapping. Its synthetic stress variants are software test cases, not experimentally failed proteins, and must not be used to report predictive accuracy.

### Protein foundation-model R&D

In August 2026, PRIB-KI completed a parameter-efficient adaptation of ESM-2 3B using experimental protein-stability labels. On the implemented held-out split with 30,000 evaluation examples, the retained result achieved:

| Target | Metric | Result |
|---|---|---:|
| Absolute stability | Spearman | **0.8875** |
| Absolute stability | Pearson | **0.8597** |
| Perturbation-induced stability change | Spearman | **0.8276** |
| Perturbation-induced stability change | Pearson | **0.8237** |
| Direction of stability change | Sign accuracy | **82.53%** |
| Embedding displacement vs. magnitude of stability change | Spearman | **0.5136** |

These are strong held-out results under the implemented split. They establish that the current research model learned substantial stability and perturbation-response signal. They do **not** yet establish generalization to unseen protein families, industrial conditions, or prospective partner panels.

The next validation gates therefore focus on leakage-aware splitting, baseline comparison, applicability domain, and prospective testing rather than simply increasing model size.

See [`research_progress_2026.md`](research_progress_2026.md) for the public technical summary.

## 2. Software validation

Every retained result should be traceable to:

- input dataset and inclusion/exclusion rules;
- source-code revision;
- dependency environment;
- model/scoring configuration;
- random seed;
- model or rule-set version;
- evaluation split definition.

The public repository should maintain deterministic output for fixed inputs and configuration, explicit invalid-sequence/missing-value handling, reproducible candidate ranking, and no silent fallback when a feature cannot be calculated.

## 3. Retrospective benchmark standard

### Primary question

For an endpoint-specific reliability model:

> Does the PRIB-KI model improve ranking of measured protein reliability or failure-related outcomes beyond simple baselines under a leakage-aware holdout?

### Candidate endpoints

Benchmarks should remain endpoint-specific. Suitable endpoints include:

- folding or thermodynamic stability;
- perturbation-induced stability change;
- expression yield;
- SEC monomer percentage;
- aggregation propensity;
- self-interaction;
- viscosity;
- solubility;
- recovery after accelerated stress.

Composite labels should only be introduced after component endpoints are understood separately.

### Dataset requirements

A benchmark dataset should include:

- traceable protein identity or sequence;
- defined assay endpoint and unit;
- assay conditions where available;
- sufficient sample size;
- documented inclusion/exclusion rules;
- compatible usage rights.

### Split strategy

Random row-level splitting is not sufficient when close sequence relatives or shared scaffolds can appear across splits. Depending on the dataset, use one or more of:

- sequence-cluster split;
- protein-family split;
- scaffold or domain split;
- temporal split;
- external holdout set.

The current August 2026 held-out result must therefore be treated as an important research milestone, not the final generalization claim.

### Baselines

At minimum compare against:

- sequence length;
- pI;
- GRAVY;
- net charge;
- hydrophobicity proxies;
- an unweighted linear score;
- a regularized linear model.

More complex models should be retained only when they improve performance under the same split and evaluation procedure.

### Metrics

For continuous endpoints:

- Spearman correlation;
- Pearson correlation where useful;
- MAE / RMSE where scale is meaningful;
- bootstrap confidence intervals.

For directional or binary endpoints:

- sign / balanced accuracy;
- ROC-AUC;
- PR-AUC;
- sensitivity and specificity.

For prioritization:

- recall of poor candidates at a fixed screening fraction;
- top-k or highest-risk-decile enrichment;
- rank correlation.

## 4. Foundation-model landscape validation

Protein-language-model work should additionally test:

- improvement over descriptor-only baselines;
- sequence/family leakage sensitivity;
- stability of learned physical-risk directions across folds;
- nearest-neighbor or embedding-distance applicability domain;
- out-of-domain flagging rather than overconfident extrapolation;
- whether embedding geometry adds information beyond the prediction head itself.

The August result currently shows a stronger learned signal for stability **direction** than for exact embedding-displacement magnitude. That distinction should remain visible in future reporting.

## 5. Prospective validation

A prospective study should follow this sequence:

1. select a new candidate panel;
2. freeze code revision, model, thresholds, and configuration;
3. generate predictions before wet-lab results are available;
4. assign candidates to predefined groups or ranking positions;
5. run predefined assays;
6. compare frozen predictions with results;
7. document incorrect predictions and protocol deviations;
8. only then update the model.

A practical first study would use **30–50 candidates** and at least two orthogonal readouts, for example expression yield plus SEC monomer percentage or DSF/nanoDSF.

Prospective validation should report top-risk enrichment, false-negative failures, calibration drift, and whether the recommended next assay would have reduced uncertainty.

## 6. TargetTrack retrospective controls

TargetTrack supports historical progression priors, data-engineering validation, and failure-taxonomy development, but it is not a direct performance-validation dataset for the public PRIB-KI scoring model.

Analyses must retain observed versus inferred stages, report trial and target levels separately, flag contradictions, and distinguish technical, nontechnical, unknown, and censored outcomes.

The current local TargetTrack snapshot contains terminal `work stopped` records without an explicit technical-failure category in the parsed status history. Those records must remain unknown rather than being silently converted into molecular failures.

## 7. Sensitivity and applicability domain

Each scoring release should include, where appropriate:

- feature/model ablation;
- threshold sensitivity;
- alternative normalization or calibration checks;
- protein-class-specific performance;
- missing-feature robustness;
- comparison of sequence-only and structure-augmented variants;
- similarity or embedding-distance based applicability-domain reporting.

Out-of-domain candidates should be flagged for review rather than assigned an overconfident score.

## 8. Claim ladder

PRIB-KI uses the following claim ladder:

1. **Workflow execution:** the pipeline runs reproducibly.
2. **Held-out research signal:** the model predicts/ranks an endpoint on a defined holdout.
3. **Leakage-aware generalization:** performance persists under sequence-family/scaffold/external holdouts.
4. **Prospective evidence:** frozen predictions are confirmed on newly measured candidates.
5. **Product evidence:** the workflow measurably improves experimental prioritization, cost, or cycle time in a partner setting.

As of September 2026, the foundation-model stability work has reached **Level 2**, while the public sequence-only demonstrator remains primarily Level 1. The active R&D program is now focused on moving the strongest model evidence toward Levels 3 and 4.
