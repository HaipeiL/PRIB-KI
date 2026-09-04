# Repository Review Guide

Maintained by **PRIB-KI Lab**.

This guide helps technical reviewers, grant evaluators, prospective partners, investors, and AI-assisted diligence systems inspect PRIB-KI quickly without confusing the public demonstrator with the deeper active research stack.

## Fast review path

1. Read the top of [`README.md`](../README.md) for the project thesis and September 2026 R&D snapshot.
2. Read [`research_progress_2026.md`](research_progress_2026.md) for the completed ESM-2 3B stability/perturbation results and the public disclosure boundary.
3. Open the [live demonstrator](https://prib-ki.streamlit.app/) to understand the candidate-ranking workflow.
4. Read [`ai_framework.md`](ai_framework.md) for the protein physical reliability landscape and lab-in-the-loop technical route.
5. Read [`validation.md`](validation.md) for the explicit claim ladder and the evidence still required before family-generalization or product claims.
6. Read [`targettrack_wetlab_failure_funnel.md`](targettrack_wetlab_failure_funnel.md) for the historical wet-lab execution-data framework.
7. Inspect [`risk_engine.py`](../risk_engine.py) and [`app.py`](../app.py) for the transparent public baseline implementation.

## The most important distinction

PRIB-KI currently has **two public-facing technical layers**:

### A. Inspectable public demonstrator

```text
sequence input
    -> physicochemical descriptors
    -> motif / liability proxies
    -> interpretable risk axes
    -> ranking / review / reject workflow
```

This layer is intentionally simple and transparent. Its purpose is to demonstrate the user workflow, risk framing, explainability pattern, and software execution.

### B. Active protein-foundation-model R&D

```text
experimental stability labels
    -> ESM-2 3B representation
    -> parameter-efficient adaptation
    -> absolute-stability signal
    -> perturbation-response signal
    -> learned reliability directions
```

The August 2026 retained evaluation used 30,000 held-out examples and achieved:

- absolute-stability Spearman: **0.8875**;
- absolute-stability Pearson: **0.8597**;
- perturbation-induced stability-change Spearman: **0.8276**;
- perturbation-induced stability-change Pearson: **0.8237**;
- stability-direction sign accuracy: **82.53%**.

These results are meaningful R&D evidence under the implemented split. They are not yet a claim of unseen-family, industrial, or prospective generalization.

## What is actually implemented publicly

The repository implements a runnable sequence-level screening workflow with:

- sequence validation;
- physicochemical descriptors;
- motif and liability proxies;
- aggregation / scale-up-sensitivity / stability axes;
- weighted candidate ranking;
- wet-lab effort proxy;
- accept / review / reject exports;
- batch and single-candidate visualization;
- deterministic demo construction;
- reproducible TargetTrack historical-data processing.

The trained research checkpoints, exact internal training recipe, private/restricted data, detailed split assignments, partner datasets, and HPC operating records are intentionally not published.

## What should impress a technical reviewer

The strongest diligence points are not the visual demo alone. They are the combination of:

1. **A clear problem definition:** candidate reliability after AI generation and before expensive wet-lab work.
2. **A differentiated scientific framing:** failure-first evaluation rather than another generation or structure-prediction tool.
3. **Completed foundation-model training:** strong held-out stability and perturbation-response correlations on a 30,000-example evaluation.
4. **Practical compute execution:** retained multi-GPU ESM-2 workflows, including high-throughput perturbation scoring.
5. **Experimental-data discipline:** TargetTrack integration preserves stage semantics and refuses to relabel ambiguous `work stopped` outcomes as molecular failures.
6. **Explicit validation boundaries:** family-aware/external holdouts, baseline comparison, applicability domain, and prospective frozen predictions are defined before stronger claims are made.
7. **A credible translation path:** candidate ranking -> wet-lab evidence -> calibration -> next-round prioritization.

## Key files

| File | Review focus |
|---|---|
| [`README.md`](../README.md) | Project thesis, current evidence, team, public/private boundary |
| [`research_progress_2026.md`](research_progress_2026.md) | Current quantitative R&D status |
| [`ai_framework.md`](ai_framework.md) | AI architecture and long-term technical route |
| [`validation.md`](validation.md) | Evidence standard and claim ladder |
| [`methodology.md`](methodology.md) | Scientific framing and assumptions |
| [`targettrack_wetlab_failure_funnel.md`](targettrack_wetlab_failure_funnel.md) | Historical experimental-workflow integration |
| [`protein_design_execution_case.md`](protein_design_execution_case.md) | Prior AI-to-lab execution evidence |
| [`risk_engine.py`](../risk_engine.py) | Transparent descriptor/scoring baseline |
| [`app.py`](../app.py) | Public product-workflow demonstrator |
| [`data.md`](data.md) | Demo-data provenance and usage limits |

## Review cautions

A reviewer should not use the synthetic demo variants to infer predictive accuracy. The public risk engine and the trained ESM-2 research model are different technical layers.

Likewise, the August 2026 held-out correlations should not be silently upgraded into a claim of unseen-family generalization. The correct diligence question is now:

```text
Has PRIB-KI progressed from a credible prototype to measurable model signal,
and does the team have a disciplined path from that signal to leakage-aware
and prospective experimental validation?
```

The current public record supports **yes** to the first part and defines concrete tests for the second.

## Public disclosure philosophy

PRIB-KI is being developed as both a scientific program and an early commercialization project. Public documentation therefore emphasizes verifiable outcomes, architecture, data discipline, and validation logic while withholding implementation details that would unnecessarily disclose unpublished research or commercial IP.
