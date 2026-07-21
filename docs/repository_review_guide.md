# Repository Review Guide

Maintained by **PRIB-KI Lab**.

This guide helps technical reviewers, grant evaluators, and early partners inspect the repository quickly.

## Fast path

1. Read the top of [`README.md`](../README.md) to understand the product scope.
2. Open the [live demonstrator](https://prib-ki.streamlit.app/) and run the built-in demo dataset.
3. Read [`long_term_vision.md`](long_term_vision.md) for the company-level two-route vision.
4. Read [`protein_design_execution_case.md`](protein_design_execution_case.md) for prior AI-to-lab execution evidence.
5. Read [`ai_framework.md`](ai_framework.md) to understand the technical AI route: protein physical reliability landscape and lab-in-the-loop physical AI.
6. Inspect [`risk_engine.py`](../risk_engine.py) for feature extraction and risk aggregation.
7. Read [`methodology.md`](methodology.md) to see the scientific assumptions.
8. Read [`validation.md`](validation.md) to see the required evidence before performance claims.

## What is implemented

The repository currently implements a sequence-only screening workflow:

```text
CSV sequence input
    -> validation
    -> physicochemical descriptors
    -> motif and liability proxies
    -> sub-risk scores
    -> Aggregation / ScaleUpSensitivity / Stability axes
    -> overall score and effort proxy
    -> accept / review / reject exports
```

The public prototype is suitable for:

- workflow demonstration;
- early hypothesis exploration;
- discussion with potential users;
- planning benchmark datasets;
- planning prospective validation.

It is not suitable for:

- clinical decisions;
- regulatory decisions;
- manufacturing release decisions;
- absolute failure-probability claims;
- reporting accuracy from the synthetic demo variants.

## How to run locally

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the app:

```bash
streamlit run app.py
```

Expected behavior:

- the app loads `data/demo_100.csv` by default;
- clicking `Evaluate` computes the risk table;
- accepted and rejected candidates appear in exportable groups;
- selecting a candidate shows risk metrics, radar plot, and sequence export.

## Key files

| File | Review focus |
|---|---|
| [`app.py`](../app.py) | User workflow, visual outputs, export behavior |
| [`risk_engine.py`](../risk_engine.py) | Descriptor calculation, pseudo-target mapping, scoring rules |
| [`build_demo_set.py`](../build_demo_set.py) | Demo dataset construction and synthetic variants |
| [`docs/long_term_vision.md`](long_term_vision.md) | Company-level vision and two-route strategy |
| [`docs/protein_design_execution_case.md`](protein_design_execution_case.md) | Prior RBX1 design competition, wet-lab testing, and post-generation reliability evidence |
| [`docs/ai_framework.md`](ai_framework.md) | AI technical stack, physical reliability landscape, and lab-in-the-loop route |
| [`docs/methodology.md`](methodology.md) | Scientific framing and implementation boundaries |
| [`docs/validation.md`](validation.md) | Benchmark, AI landscape validation, and prospective validation requirements |
| [`docs/data.md`](data.md) | Data provenance and usage limits |
| [`docs/roadmap.md`](roadmap.md) | Technical risk-reduction milestones |

## Review cautions

The current scoring layer uses heuristic pseudo-targets for demonstrator behavior. This is acceptable for showing a decision workflow, but not for making predictive-performance claims.

The central diligence question is therefore not "is this already a validated predictor?" but:

```text
Is the workflow credible, inspectable, and ready to generate validation evidence?
```

The planned validation work should answer whether high-risk PRIB-KI rankings enrich experimentally problematic candidates under predefined assay endpoints.

For the longer AI route, reviewers should additionally ask whether protein language-model embeddings and failure-anchor rules improve over simple descriptor baselines, whether out-of-domain candidates are flagged, and whether frozen wet-lab feedback can update the risk landscape reproducibly.

The prior AI protein-design sprint should be reviewed as team execution evidence. It should not be counted as validation evidence for the current PRIB-KI scoring model unless the exact design records, wet-lab readouts, and frozen prediction protocol are separately documented.
