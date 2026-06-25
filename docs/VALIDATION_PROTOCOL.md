# PRIB-KI Validation Protocol

## 1. Objective

This protocol defines how PRIB-KI should progress from a research demonstrator to an experimentally supported decision tool.

The protocol separates engineering reproducibility, retrospective benchmarking, and prospective wet-lab validation. Performance claims should only be made at the level supported by completed validation.

## 2. Validation Principles

- Predictions must be frozen before evaluation.
- Training, calibration, and test data must be separated.
- Closely related sequences must not be split across train and test sets without explicit justification.
- Synthetic stress variants may be used for software and sensitivity testing, but not as experimental ground truth.
- Absolute risk claims require calibration against experimental outcomes.
- Uncertainty and applicability domain must be reported together with scores.
- Negative results and failure cases must be documented.

## 3. Level 1 — Engineering and Reproducibility Validation

### 3.1 Goals

Confirm that the software behaves deterministically, handles inputs correctly, and produces traceable outputs.

### 3.2 Required checks

- identical input and configuration produce identical output;
- random seeds are fixed and recorded;
- sequence validation rejects malformed inputs;
- descriptor calculations pass unit tests;
- missing values and edge cases are handled explicitly;
- score ranges and thresholds are tested;
- configuration changes are logged;
- data and model versions are included in exported results.

### 3.3 Acceptance criteria

- all automated tests pass;
- no undocumented nondeterminism;
- every output can be linked to an input file, software version, and configuration;
- repeated execution produces the same ranking within numerical tolerance.

## 4. Level 2 — Retrospective Public-Data Validation

### 4.1 Primary hypothesis

Higher PRIB-KI risk scores are associated with poorer experimentally measured developability outcomes.

### 4.2 Candidate endpoints

- expression yield;
- SEC monomer percentage;
- aggregation propensity;
- thermal stability;
- self-interaction;
- viscosity;
- accelerated degradation;
- solubility;
- recovery after stress.

Each endpoint must be analyzed separately before any composite label is introduced.

### 4.3 Dataset inclusion criteria

A dataset may be included when it provides:

- traceable protein identity or sequence;
- clearly defined experimental endpoint;
- documented assay conditions;
- sufficient sample size for the intended analysis;
- licensing or usage conditions compatible with the project.

### 4.4 Exclusion criteria

Exclude or separately flag:

- unresolved sequence identity;
- duplicated measurements without aggregation rules;
- incompatible assay conditions pooled without correction;
- ambiguous labels;
- samples with excessive missingness;
- data used during model development when evaluating final performance.

### 4.5 Leakage control

Use one or more of the following, depending on dataset size:

- cluster split by sequence identity;
- protein-family split;
- scaffold or domain split;
- temporal split;
- external dataset holdout.

Random row-level splitting alone is not sufficient where homologous proteins are present.

### 4.6 Baselines

PRIB-KI should be compared against simple and transparent baselines:

- sequence length;
- pI;
- GRAVY or mean hydropathy;
- net charge;
- hydrophobic patch proxy;
- unweighted linear score;
- standard regularized regression or classification model.

A more complex model should only be retained when it provides reproducible improvement over these baselines.

### 4.7 Evaluation metrics

For continuous endpoints:

- Spearman correlation;
- Pearson correlation where justified;
- mean absolute error;
- root mean squared error;
- confidence intervals from bootstrap resampling.

For binary or ordinal endpoints:

- ROC-AUC;
- PR-AUC;
- balanced accuracy;
- sensitivity and specificity;
- top-k enrichment;
- calibration curve and Brier score when probability outputs are used.

For ranking use cases:

- rank correlation;
- enrichment of experimentally poor candidates in the highest-risk fraction;
- recall of poor candidates at fixed screening capacity.

### 4.8 Predefined success criteria

Success criteria must be selected before final evaluation. An example for an early prototype is:

- statistically supported association with at least one independent experimental endpoint;
- improvement over simple physicochemical baselines;
- stable performance across bootstrap resamples;
- documented failure cases;
- no major collapse under protein-family-aware splitting.

These criteria are placeholders until fixed for a specific benchmark study.

## 5. Level 3 — Prospective Blinded Validation

### 5.1 Study design

1. Select a new candidate panel not used in model development.
2. Record the intended protein class and use context.
3. Freeze model version, configuration, and thresholds.
4. Generate predictions before wet-lab results are available.
5. Assign candidates to low-, intermediate-, and high-risk groups.
6. Perform predefined assays using standardized protocols.
7. Compare predictions with experimental outcomes.
8. Publish both successful and unsuccessful predictions internally or externally.

### 5.2 Suggested initial panel

- 30–50 candidates;
- sufficient diversity to avoid near-duplicate sequences;
- balanced representation across predicted risk groups;
- at least two orthogonal experimental endpoints.

### 5.3 Suggested assay package

A practical initial package may include:

- expression yield;
- SEC monomer percentage;
- nanoDSF or DSF;
- aggregation onset or turbidity;
- one self-interaction or solubility-related assay.

### 5.4 Primary prospective endpoint

The primary endpoint should be fixed before testing. A suitable early endpoint could be enrichment of experimentally poor candidates in the predefined high-risk group.

## 6. Sensitivity and Ablation Analysis

The following analyses should accompany model validation:

- removal of each risk module;
- alternative feature normalizations;
- threshold sensitivity;
- weight sensitivity;
- protein-class-specific performance;
- performance under missing descriptors;
- comparison of sequence-only and sequence-plus-structure models.

## 7. Applicability Domain

Each prediction should eventually report whether the candidate is within the model's supported domain.

Potential domain indicators include:

- sequence length range;
- protein class;
- similarity to validated data;
- descriptor-space distance;
- missing contextual information;
- unusual composition or architecture.

Out-of-domain candidates should be marked for review rather than assigned overconfident scores.

## 8. Reporting Standard

Every validation report should include:

- dataset version and source;
- inclusion and exclusion flow;
- endpoint definition;
- split strategy;
- model and configuration version;
- baseline results;
- confidence intervals;
- calibration status;
- known limitations;
- complete list of deviations from the predefined protocol.

## 9. Current Validation Status

The current repository demonstrates workflow execution and controlled descriptor stratification. It has not yet established experimentally validated predictive performance.

Synthetic stress variants are used only to test demonstrator behavior and must not be interpreted as experimentally failed proteins.