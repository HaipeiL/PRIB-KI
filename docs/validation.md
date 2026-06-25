# Validation

Maintained by **PRIB-KI Lab**.

This document defines the minimum validation work required before performance claims are made.

## 1. Software validation

The repository must pass the following checks:

- deterministic output for fixed input, seed, and configuration;
- unit tests for descriptor calculations;
- explicit handling of invalid sequences and missing values;
- reproducible candidate ranking;
- versioned configuration and output metadata;
- no silent fallback when a feature cannot be calculated.

A result should be traceable to:

- input dataset;
- source code revision;
- dependency environment;
- scoring configuration;
- random seed;
- model or rule-set version.

## 2. Retrospective benchmark

### Primary question

Does the PRIB-KI ranking enrich experimentally poor candidates at the high-risk end of the distribution?

### Candidate endpoints

Benchmarks should be endpoint-specific. Suitable endpoints include:

- expression yield;
- SEC monomer percentage;
- thermal stability;
- aggregation propensity;
- self-interaction;
- viscosity;
- solubility;
- recovery after accelerated stress.

Composite labels should only be introduced after the component endpoints are analyzed separately.

### Dataset requirements

A benchmark dataset must include:

- traceable protein identity or sequence;
- defined assay endpoint;
- assay conditions;
- sufficient sample size;
- documented inclusion and exclusion rules;
- compatible usage rights.

### Split strategy

Random row-level splitting is not sufficient when related proteins are present. Depending on the dataset, use:

- sequence-cluster split;
- protein-family split;
- scaffold or domain split;
- temporal split;
- external holdout set.

### Baselines

At minimum, compare against:

- sequence length;
- pI;
- GRAVY;
- net charge;
- a hydrophobicity proxy;
- an unweighted linear score;
- a regularized linear model.

A more complex model should only be retained if it improves performance under the same split and evaluation procedure.

### Metrics

For continuous endpoints:

- Spearman correlation;
- mean absolute error;
- bootstrap confidence intervals.

For binary or ordinal endpoints:

- ROC-AUC;
- PR-AUC;
- balanced accuracy;
- sensitivity and specificity;
- top-k enrichment.

For prioritization:

- recall of experimentally poor candidates at a fixed screening fraction;
- enrichment in the highest-risk decile or quartile;
- rank correlation.

## 3. Prospective validation

A prospective study should follow this sequence:

1. select a new candidate panel;
2. freeze the code revision, model, and thresholds;
3. generate predictions before wet-lab results are available;
4. assign candidates to predefined risk groups;
5. run predefined assays;
6. compare predictions with results;
7. record failed predictions and protocol deviations.

A practical first study would use 30–50 candidates and at least two orthogonal readouts, for example:

- expression yield;
- SEC monomer percentage;
- nanoDSF or DSF;
- aggregation onset or turbidity;
- one self-interaction or solubility assay.

## 4. Sensitivity analysis

The following checks are required for each scoring release:

- feature ablation;
- weight sensitivity;
- threshold sensitivity;
- alternative normalization methods;
- protein-class-specific performance;
- missing-feature robustness;
- comparison of sequence-only and structure-augmented variants.

## 5. Applicability domain

Predictions should eventually report whether a candidate lies inside the validated domain. Relevant factors include:

- sequence length;
- protein class;
- similarity to benchmark data;
- descriptor-space distance;
- missing context;
- unusual composition or architecture.

Out-of-domain candidates should be flagged for review rather than assigned an overconfident score.

## 6. Current status

The current repository demonstrates workflow execution and controlled descriptor stratification.

The synthetic stress variants in the demo dataset are software test cases. They are not experimentally failed proteins and must not be used to report predictive accuracy.