# Methodology

Maintained by **PRIB-KI Lab**.

## Scope

PRIB-KI is an early-stage screening workflow for protein developability and reliability assessment. The current implementation operates primarily on sequence-derived descriptors and produces relative risk rankings for candidate prioritization.

The repository should be read as a research prototype. It is not a validated predictor of clinical or manufacturing failure and is not intended to replace experimental developability assays.

## Working model

The pipeline follows a failure-oriented workflow:

```text
sequence input
    -> input validation
    -> descriptor calculation
    -> failure-mode features
    -> risk aggregation
    -> candidate ranking
    -> review / retain / deprioritize
```

The central assumption is that a structurally plausible protein can still carry condition-dependent liabilities that become relevant during expression, purification, formulation, storage, concentration, transport, or repeated handling.

## AI framework

PRIB-KI uses AI as a physical reliability layer rather than as an autonomous protein generator. The intended method has two stages.

Stage A builds a **protein physical reliability landscape**:

```text
amino-acid sequence
    -> protein language-model embedding
    -> optional structure and condition augmentation
    -> learned physical-risk directions
    -> danger-zone and applicability-domain rules
    -> candidate reliability ranking
```

Stage B builds a **lab-in-the-loop physical AI autopilot**:

```text
current risk landscape
    -> next experiment selection
    -> frozen wet-lab prediction
    -> assay readout ingestion
    -> calibration update
    -> next-round reliability map
```

This makes the project part of three overlapping categories:

- **Drug x AI:** a candidate-prioritization layer between computational design and developability experiments.
- **Physical AI:** a hybrid model that keeps biophysical failure mechanisms visible inside the scoring logic.
- **AI for life science:** a design-test-learn workflow where model predictions are updated only after traceable experimental feedback.

The full architecture and technical route are described in [`ai_framework.md`](ai_framework.md).

## AI methods used in the target system

| Method | Role in PRIB-KI |
|---|---|
| Protein language-model embedding | Maps amino-acid sequences into an AI-learned protein landscape |
| Structure and condition augmentation | Adds geometry and stress context to the candidate representation |
| Risk-direction learning | Learns aggregation, stability, expression, SEC, solubility, or self-interaction directions from endpoint labels |
| Failure-anchor retrieval | Compares new candidates against experimentally stable or fragile anchors |
| Applicability-domain detection | Identifies candidates outside validated evidence rather than forcing a confident score |
| Active learning and Bayesian experimental design | Selects the next candidate panel or stress condition that should most improve the model |
| Model calibration | Updates risk estimates after frozen predictions are compared with wet-lab readouts |

## Reliability dimensions

The project currently separates the following concepts:

- thermodynamic stability;
- kinetic stability;
- colloidal stability;
- chemical stability;
- conformational robustness;
- process sensitivity;
- mechanical or pathway-dependent robustness.

These dimensions are related but not interchangeable. A proxy used for one dimension must not be presented as direct evidence for another.

## Evidence mapping

Each feature should be traceable through the following chain:

```text
perturbation
    -> molecular response
    -> failure mechanism
    -> experimental readout
    -> computational proxy
    -> decision signal
```

Current examples:

| Condition | Failure mode | Typical experimental readout | Current implementation |
|---|---|---|---|
| elevated temperature | partial unfolding or aggregation | Tm, Tagg, SEC monomer loss | sequence-level proxy |
| agitation or interfaces | interface-induced aggregation | turbidity, particles, monomer loss | hydrophobicity and charge proxies |
| high concentration | self-association or viscosity | kD, AC-SINS, viscosity | charge and hydrophobicity proxies |
| oxidative stress | Met/Trp oxidation | LC-MS degradation profile | residue-level proxy |
| pH shift | solubility or colloidal instability | recovery, aggregation, solubility | pI and charge proxies |
| mechanical loading | geometry-dependent unfolding or rupture | AFM-SMFS, force-clamp, SMD | not implemented |

## Current implementation

The current codebase provides:

- sequence validation;
- physicochemical descriptor calculation;
- motif-based liability features;
- relative risk-axis calculation;
- AI-weighted pseudo-target mapping for demonstrator behavior;
- configurable score aggregation;
- batch ranking and visualization.

The current score is a relative screening index. It must not be interpreted as an experimentally calibrated probability.

## Historical progression priors

The TargetTrack data module adds a separate, descriptive retrospective layer
for historical structural-genomics workflows. It represents a progression
history as events, separates direct observations from inferred prerequisites,
and keeps technical failures distinct from nontechnical and unknown project
stops. Trial/construct units and target units are reported separately so a
failed trial is not erased by a later successful retry.

This is deliberately not connected to the demonstrator score. Before any
TargetTrack-derived supervised model is trained, raw mappings must be reviewed,
the descriptive funnel must pass audit, and labels must exclude unknown,
nontechnical and contradictory outcomes. Details are in
[targettrack_wetlab_failure_funnel.md](targettrack_wetlab_failure_funnel.md).

The current local TargetTrack snapshot run has a complete mapped status
inventory, but no explicit technical terminal-failure labels in the parsed
history. Its \`work stopped\` records remain unknown outcomes, so the source
supports descriptive progression and taxonomy work only—not a supervised
failure model or resolved-outcome probability estimate.

## Planned extensions

The next technical steps are:

- fixed feature and configuration schemas;
- benchmark data schemas with assay condition metadata;
- protein language-model embeddings where benchmark evidence supports them;
- structure-derived descriptors and confidence metrics where justified;
- condition-specific inputs for pH, temperature, concentration, buffer, storage, and stress context;
- stable and failed anchors in embedding space;
- learned physical-risk directions and danger-zone rules;
- experimentally annotated benchmark datasets;
- uncertainty estimates;
- applicability-domain detection;
- active-learning logic for selecting the next wet-lab panel;
- frozen prediction and evidence registry;
- blinded prospective validation.

## Design rules

1. Keep feature extraction separate from scoring.
2. Store weights and thresholds in versioned configuration files.
3. Record dataset, software, and configuration versions with every result.
4. Do not use synthetic variants as experimental ground truth.
5. Compare new models against simple physicochemical baselines.
6. Report out-of-domain cases instead of forcing a confident score.
7. Treat wet-lab validation as the final source of evidence.
