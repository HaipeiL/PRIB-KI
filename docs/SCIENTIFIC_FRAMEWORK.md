# PRIB-KI Scientific Framework

## 1. Purpose

PRIB-KI is a research-stage framework for early, failure-oriented assessment of protein developability and industrial reliability.

The central question is not only whether a protein can adopt a plausible structure, but whether it is likely to remain usable under non-ideal conditions relevant to development, production, formulation, storage, transport, or repeated handling.

PRIB-KI is intended as a decision-support layer upstream of costly experimental validation. It does not replace wet-lab assays and does not currently provide clinically validated failure probabilities.

## 2. Core Hypothesis

Structurally plausible proteins may still fail because of latent, condition-dependent failure pathways that are not captured by structure plausibility alone.

PRIB-KI therefore follows a failure-first logic:

1. define relevant perturbations;
2. identify plausible molecular failure mechanisms;
3. connect those mechanisms to measurable assay endpoints;
4. derive computational proxies;
5. integrate the evidence into an interpretable risk profile;
6. support retain, review, or deprioritize decisions.

## 3. Reliability Dimensions

Protein reliability is treated as a multidimensional concept.

### 3.1 Thermodynamic stability
The equilibrium preference for folded versus unfolded states.

### 3.2 Kinetic stability
The height of barriers that delay unfolding, dissociation, or irreversible degradation.

### 3.3 Colloidal stability
The tendency to remain dispersed rather than self-associate, aggregate, or phase-separate.

### 3.4 Chemical stability
Sensitivity to oxidation, deamidation, isomerization, clipping, hydrolysis, or other covalent degradation processes.

### 3.5 Conformational robustness
Resistance to excessive flexibility, partial unfolding, or stress-induced population shifts.

### 3.6 Mechanical and pathway-dependent robustness
Sensitivity to force geometry, loading history, interfaces, flow, agitation, or other non-equilibrium perturbations.

### 3.7 Process robustness
Tolerance to expression, purification, concentration, formulation, storage, and scale-up conditions.

## 4. Scientific Causal Chain

```text
Perturbation
    ↓
Molecular response
    ↓
Failure mechanism
    ↓
Experimental endpoint
    ↓
Computational proxy
    ↓
Risk evidence
    ↓
Decision recommendation
```

This causal chain is used to prevent descriptors from being treated as direct proof of failure. Each proxy must be linked to a stated hypothesis, an experimental endpoint, and a validation plan.

## 5. Current Failure-Mode Map

| Perturbation | Potential failure mechanism | Experimental endpoint | Current computational proxy status |
|---|---|---|---|
| Elevated temperature | partial unfolding, aggregation | Tm, Tagg, onset temperature, SEC loss | sequence-level proxy only |
| Agitation or interfaces | interface-induced aggregation | particle count, turbidity, monomer loss | hydrophobicity and charge proxies |
| Concentration | self-association, viscosity | kD, AC-SINS, viscosity | charge and hydrophobicity proxies |
| Oxidative stress | Met/Trp oxidation | LC-MS degradation profile | residue and motif proxy |
| pH shift | colloidal instability | solubility, aggregation, recovery | pI and charge proxy |
| Repeated handling | pathway-dependent instability | accelerated stress assays | not yet explicitly modeled |
| Mechanical loading | geometry-dependent unfolding or rupture | AFM-SMFS, force-clamp, SMD | future module |

## 6. Representation Strategy

PRIB-KI is designed to combine several representation levels:

- sequence-derived descriptors;
- structure-derived descriptors;
- perturbation metadata;
- application and formulation context;
- experimentally observed developability outcomes.

The current demonstrator is primarily sequence-based. Structural, perturbation-specific, and experimentally calibrated representations are part of the development roadmap.

## 7. Evidence Hierarchy

PRIB-KI distinguishes between four levels of evidence:

1. **Mechanistic evidence** — supported by established biophysical principles.
2. **Proxy evidence** — indirect computational indicators linked to a failure hypothesis.
3. **Retrospective evidence** — association with annotated experimental datasets.
4. **Prospective evidence** — blinded predictions tested in new experiments.

Only the last two levels can support quantitative performance claims.

## 8. Current Scope

The current version demonstrates:

- reproducible sequence validation;
- interpretable descriptor extraction;
- heuristic and model-assisted risk aggregation;
- candidate ranking;
- visual decision support.

## 9. Current Limitations

The current version does not yet:

- estimate absolute probabilities of development failure;
- model all protein classes equally well;
- incorporate formulation-specific conditions in a calibrated manner;
- replace experimental developability assays;
- support clinical, regulatory, or manufacturing release decisions;
- demonstrate prospective predictive validity.

## 10. Research Direction

The long-term objective is to develop a perturbation-aware reliability layer that combines biophysical constraints, experimental evidence, uncertainty estimation, and application-specific context.

The intended output is not a single universal score, but an auditable risk profile that identifies likely failure drivers and recommends the most informative next validation step.