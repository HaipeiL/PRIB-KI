# Long-term Vision

Maintained by **PRIB-KI Lab**.

This document states the long-term company vision for PRIB-KI. It is intentionally more strategic than the implementation docs. Current product claims remain bounded by the validation plan in [`validation.md`](validation.md).

## Core belief

PRIB-KI is built on a 10-year judgement: AI protein engineering is moving from a **generation bottleneck** to a **reliability bottleneck**.

The first wave of AI-biology systems increases the number of plausible protein candidates that can be generated, folded, scored, or proposed. The next industrial bottleneck is different: deciding which candidates deserve scarce wet-lab attention, formulation work, scale-up learning, and scientific review.

In that world, the strategic asset is not only a generator. It is a system that can answer:

```text
Which AI-designed proteins are physically reliable enough to earn the next experiment?
```

PRIB-KI aims to become that physical reliability layer.

## Execution basis

The vision is informed by a prior AI-to-lab execution case documented in [`protein_design_execution_case.md`](protein_design_execution_case.md). In the GEM x Adaptyv RBX1 Binder Design Competition linked to the ICLR 2026 GEM Workshop, the team used modern protein-design tools, including Proteina-Complexa and BindCraft, applied reliability-oriented candidate triage, and advanced 7 of 21 submitted candidates to wet-lab testing. This corresponded to approximately the top 2.5% of more than 12,000 global submissions by external selection performance.

The campaign did not ultimately produce a confirmed functional RBX1 binder. That result matters because it exposed the practical bottleneck PRIB-KI is built around: after AI tools generate plausible proteins, the decisive question becomes which candidates deserve experimental resources.

This case is not treated as PRIB-KI model validation. It is execution evidence that the AI-to-lab workflow is not hypothetical for the project.

## The two routes

The long-term vision has two connected AI routes.

### Route A - Build the map

**Protein physical reliability landscape**

Route A turns amino-acid sequences into an AI-learned protein landscape and annotates that landscape with physical danger zones.

```text
AA sequence
    -> protein language-model embedding
    -> structure and condition augmentation
    -> physical-risk directions
    -> danger-zone rules
    -> reliability ranking
```

Technical stack:

| Layer | Methods | Purpose |
|---|---|---|
| Protein embedding | ESM-family, ProtT5-family, or comparable protein language-model embeddings | Place candidates in a learned protein landscape |
| Physical augmentation | Predicted structure, confidence metrics, surface exposure, charge/hydrophobic patches, condition metadata | Add geometry and stress context |
| Risk-direction learning | Linear probes, regularized models, gradient-boosted trees, metric learning, contrastive learning | Learn directions linked to aggregation, instability, poor expression, poor SEC, solubility, or self-interaction |
| Danger-zone rules | Failure-anchor retrieval, clustering, calibrated thresholds, applicability-domain checks | Convert representation space into reviewable reliability decisions |

The key insight is that embedding space alone is not the asset. The asset is an AI protein landscape with experimentally grounded physical-risk directions.

Route A should be proven with:

- improvement over physicochemical baselines such as length, pI, charge, GRAVY, and hydrophobicity;
- top-k enrichment of experimentally poor candidates;
- separation between stable anchors and fragile anchors under leakage-aware splits;
- calibrated danger-zone thresholds;
- explicit out-of-domain flags.

### Route B - Navigate the map

**Lab-in-the-loop physical AI autopilot**

Route B turns the reliability landscape into a closed learning system. It decides what to test next, records frozen predictions, learns from wet-lab outcomes, and updates the map.

```text
risk landscape
    -> next experiment policy
    -> wet-lab / robot / partner assay
    -> frozen evidence registry
    -> calibration update
    -> next-round selection
```

Technical stack:

| Layer | Methods | Purpose |
|---|---|---|
| Experiment selection | Active learning, uncertainty sampling, diversity-aware batch selection | Choose candidate panels that improve the model fastest |
| Condition navigation | Bayesian experimental design or Bayesian optimization | Explore pH, temperature, concentration, buffer, storage, and stress conditions |
| Evidence update | Calibration models, conformal prediction, uncertainty updates, error analysis | Turn assay readouts into better risk estimates |
| AI agents | Assay-planning agent, evidence-tracking agent, report-generation agent | Operate the scientific workflow without replacing human review |
| Automation layer | ELN/LIMS integration, robot-ready assay manifests, API workflows | Connect prediction, experiment execution, and evidence capture |

The key insight is that agents should not be framed as autonomous scientists. They should be evidence-loop operators that keep predictions, assays, data, and model updates aligned.

Route B should be proven with:

- prospective enrichment of predefined wet-lab failures;
- reduced low-value screening effort in partner panels;
- shorter candidate-to-report turnaround;
- traceability from every output to data, model version, configuration, and assay context;
- reproducible calibration updates after frozen predictions are compared with readouts.

## Why this can matter

If generative protein AI makes plausible candidate creation abundant, then the bottleneck moves downstream. Teams will need a reliability infrastructure that sits between design and experiment:

```text
generation abundance
    -> selection pressure
    -> physical reliability intelligence
    -> lab-in-the-loop learning
```

PRIB-KI's long-term goal is to change protein engineering decisions from:

```text
Which design looks plausible?
```

to:

```text
Which design is physically reliable enough to spend the next experiment on?
```

That is the project's potential non-replaceable value in AI protein design and drug x AI workflows.

## Time horizon

| Horizon | Strategic objective | Evidence required |
|---|---|---|
| 0-2 years | Build the transparent reliability workflow and first benchmark-ready AI landscape | Reproducible demo, endpoint datasets, baselines, embedding-landscape validation |
| 3-5 years | Prove physical-risk maps on partner panels | Prospective wet-lab validation, calibrated risk directions, out-of-domain detection, assay recommendation |
| 5-10 years | Operate as a lab-in-the-loop reliability infrastructure | Evidence registry, active learning, automation connectors, repeatable partner workflows, compounding validation data |

## Boundaries

This vision is not a current performance claim. The public prototype is a research demonstrator with sequence-derived descriptors and relative risk rankings. The long-term vision becomes credible only if each stage passes the evidence gates described in [`validation.md`](validation.md) and [`roadmap.md`](roadmap.md).
