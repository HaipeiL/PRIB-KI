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
- configurable score aggregation;
- batch ranking and visualization.

The current score is a relative screening index. It must not be interpreted as an experimentally calibrated probability.

## Planned extensions

The next technical steps are:

- fixed feature and configuration schemas;
- structure-derived descriptors;
- condition-specific inputs;
- experimentally annotated benchmark datasets;
- uncertainty estimates;
- applicability-domain detection;
- blinded prospective validation.

## Design rules

1. Keep feature extraction separate from scoring.
2. Store weights and thresholds in versioned configuration files.
3. Record dataset, software, and configuration versions with every result.
4. Do not use synthetic variants as experimental ground truth.
5. Compare new models against simple physicochemical baselines.
6. Report out-of-domain cases instead of forcing a confident score.
7. Treat wet-lab validation as the final source of evidence.