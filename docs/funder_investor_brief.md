# Funder and Investor Brief

Maintained by **PRIB-KI Lab**.

This document is a concise diligence guide for government grant reviewers, university transfer offices, angel investors, and strategic partners.

## Executive snapshot

| Item | Summary |
|---|---|
| Project | PRIB-KI |
| Category | AI-enabled protein reliability and developability screening |
| Stage | Research prototype with live Streamlit demonstrator |
| Wedge | Candidate prioritization between AI design and wet-lab validation |
| Buyer pain | Too many plausible protein candidates and too little validation budget |
| Initial users | AI protein design teams, biotech/pharma R&D, enzyme-engineering groups, CRO/CDMO teams, academic translation groups |
| Current evidence | Workflow demonstrator, deterministic demo dataset, documented method boundaries |
| Next evidence step | Endpoint-specific benchmarks and prospective wet-lab validation |

## Problem

AI-enabled protein design is accelerating candidate generation. Wet-lab validation capacity, however, remains expensive, slow, and limited. Many candidates can look plausible by sequence, structure, or generation score, yet later fail because of aggregation, poor expression, instability, conformational lability, formulation sensitivity, or scale-up stress.

This creates a decision-quality bottleneck:

```text
more generated candidates
        -> limited validation capacity
        -> expensive triage decisions
        -> late physical failure if weak candidates enter experiments
```

PRIB-KI focuses on reducing the cost of choosing the wrong candidates early.

## Solution thesis

PRIB-KI is a failure-first reliability layer. Instead of asking only whether a protein can be generated or structurally predicted, the platform asks whether a candidate appears physically reliable enough to justify further experimental investment.

The prototype currently:

- ingests candidate sequences;
- extracts interpretable descriptors;
- maps descriptors to reliability-relevant risk axes;
- ranks candidates for accept, review, or reject decisions;
- exports decision tables for wet-lab planning.

The intended future product is a validated decision-support platform that integrates with existing design, assay-planning, and portfolio-review workflows.

## Why this is differentiated

PRIB-KI is not positioned as another protein design engine. Its differentiation is the question it asks:

| Conventional design workflow | PRIB-KI |
|---|---|
| Can we generate a plausible candidate? | Should we spend experiments on this candidate? |
| Focus on success metrics | Focus on failure pathways and reliability risk |
| Often structure- or likelihood-centered | Perturbation-aware and developability-centered |
| Produces candidates | Filters and prioritizes candidates |

The project is grounded in protein biophysics, single-molecule failure logic, physical modeling, and decision-oriented software workflow.

## Current technical status

Implemented in the public repository:

- Streamlit review interface;
- deterministic demo dataset;
- sequence validation;
- physicochemical descriptors;
- motif-based liability proxies;
- aggregation, scale-up sensitivity, and stability axes;
- overall risk aggregation;
- wet-lab effort index as a planning proxy;
- exportable accept/review/reject groups;
- documentation of assumptions, data limits, roadmap, and validation plan.

Current limitation:

The public scores are relative screening indices, not experimentally calibrated probabilities. The demo data includes synthetic stress variants for software testing, not experimental negative labels.

## Evidence standard

PRIB-KI should be evaluated using the following evidence ladder:

1. **Software reproducibility**  
   Fixed inputs, fixed configuration, fixed seed, deterministic rankings, and documented output provenance.

2. **Retrospective endpoint benchmarks**  
   Compare risk rankings with measured endpoints such as expression yield, SEC monomer percentage, thermal stability, aggregation onset, solubility, self-interaction, or viscosity.

3. **Baseline comparison**  
   Compare against length, pI, GRAVY, net charge, hydrophobicity proxies, and simple linear scores before claiming model advantage.

4. **Family-aware evaluation**  
   Avoid row-level leakage by using sequence-cluster, scaffold, domain, protein-family, temporal, or external holdout splits.

5. **Prospective validation**  
   Freeze code, configuration, and thresholds before wet-lab results are known; then test whether high-risk candidates are enriched for experimental problems.

## Proposed 12-month de-risking plan

| Period | Milestone | Diligence value |
|---|---|---|
| Months 1-3 | Reproducible build, tests, schemas, output provenance | Confirms engineering reliability |
| Months 3-6 | Endpoint-specific public benchmarks and baselines | Tests scientific signal beyond simple descriptors |
| Months 6-9 | Wet-lab partner and frozen prospective panel | Converts prototype into evidence-generating system |
| Months 9-12 | Prospective report, calibration, pilot package | Supports grant continuation, seed readiness, or paid pilots |

## Commercialization path

The likely commercialization sequence is:

1. **Paid validation pilot**  
   Customer or partner candidate panel, documented assay endpoints, before/after workflow analysis.

2. **Software license or SaaS workflow**  
   Batch screening, project-specific configuration, export reports, and portfolio review.

3. **API and enterprise integration**  
   Integration with design pipelines, ELN/LIMS, assay planning, and data governance.

4. **Validated specialized modules**  
   Endpoint-specific and protein-class-specific risk modules after benchmark evidence supports them.

## Defensibility

Potential defensibility should come from:

- curated benchmark and validation datasets;
- experimentally linked failure-mode modules;
- versioned scoring configurations;
- partner-specific feedback loops;
- domain expertise in protein mechanics and biophysical failure;
- practical integration into wet-lab decision workflows.

The public repository is a demonstrator and transparency layer. Proprietary partner datasets, calibrated modules, and validated workflow integrations can be managed separately.

## Key risks and mitigation

| Risk | Mitigation |
|---|---|
| Scores remain heuristic | Build endpoint-specific benchmarks and prospective validation |
| Data leakage in related protein families | Use family-aware and cluster-aware splits |
| Claims exceed evidence | Maintain explicit claims policy and validation gates |
| Adoption blocked by workflow friction | Provide CSV/API exports, audit-ready reports, and interpretable risk drivers |
| Overgeneralization across protein classes | Add applicability-domain flags and class-specific modules |
| Wet-lab validation cost | Start with a focused 30-50 candidate panel and two orthogonal readouts |

## What reviewers should look for

Reviewers can inspect:

- [`README.md`](../README.md) for the project overview;
- [`risk_engine.py`](../risk_engine.py) for implemented scoring logic;
- [`app.py`](../app.py) for the user-facing workflow;
- [`docs/methodology.md`](methodology.md) for scientific assumptions;
- [`docs/validation.md`](validation.md) for the evidence plan;
- [`docs/data.md`](data.md) for data provenance and limits;
- [`docs/roadmap.md`](roadmap.md) for technical milestones.

## Claims policy

Acceptable current claims:

- PRIB-KI is a research prototype.
- The prototype demonstrates an end-to-end protein candidate triage workflow.
- Current outputs are relative screening indices.
- The project has a documented validation plan.

Claims that should wait for validation:

- PRIB-KI predicts clinical or manufacturing failure.
- PRIB-KI reduces wet-lab cost by a fixed percentage across all protein classes.
- PRIB-KI outperforms baselines on experimental endpoints.
- PRIB-KI replaces developability assays.

## Contact

PRIB-KI Lab  
Dr. Haipei Liu  
haipei.thu@gmail.com
