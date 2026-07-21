# AI Framework and Technical Route

Maintained by **PRIB-KI Lab**.

This document explains how PRIB-KI should be understood as an AI for life science project, how it relates to drug x AI and physical AI, and how the technical route should evolve from the current demonstrator into a validated reliability platform.

For the company-level long-term vision, see [`long_term_vision.md`](long_term_vision.md).

## One-line positioning

PRIB-KI is a **physical-AI reliability layer for AI-designed proteins**: it does not generate proteins; it helps decide which generated or discovered protein candidates are physically reliable enough to justify experimental budget.

```text
AI generation / protein design
        -> PRIB-KI physical reliability screening
        -> wet-lab validation planning
        -> experimental feedback and model calibration
```

## Why this is AI for life science

Life-science AI is shifting from single-task prediction toward full design-test-learn workflows. Structure predictors, protein language models, inverse-folding models, and diffusion-based design systems can now create or evaluate large candidate libraries. This creates a new bottleneck: too many plausible candidates compete for limited expression, purification, formulation, stability, and scale-up experiments.

PRIB-KI addresses this bottleneck by adding a failure-oriented selection layer. Its AI value is not "make another protein"; its value is to combine learned biological representations with biophysical failure logic so teams can spend wet-lab capacity on better-ranked candidates.

## Three complementary framings

| Framing | What it means for PRIB-KI | What PRIB-KI should not claim yet |
|---|---|---|
| Drug x AI | Candidate triage between computational design and wet-lab developability testing for therapeutic proteins, antibodies, enzymes, and engineered binders. | Clinical candidate selection, efficacy prediction, regulatory decision support, or manufacturing release decisions. |
| Physical AI | Hybrid AI that maps learned protein representations to physical danger zones such as aggregation, instability, self-interaction, viscosity, chemical liability, formulation stress, and process sensitivity. | A black-box score with no mechanistic interpretation or no link to experimental readouts. |
| AI for life science | A design-test-learn platform that can connect sequence, structure, foundation-model embeddings, assay endpoints, and partner feedback. | A universal model that works across every protein class before benchmark evidence exists. |

## Two-system AI route

PRIB-KI's AI route is intentionally concentrated into two systems.

### A. Protein physical reliability landscape

Goal: build an AI-learned landscape where protein candidates can be located by sequence, structure, and condition, then interpreted through physical danger zones.

```text
amino-acid sequence
    -> protein language-model embedding
    -> optional structure and condition augmentation
    -> physical-risk directions
    -> danger-zone rules
    -> reliability ranking
```

Key technical stack:

| Component | AI method | Technical role |
|---|---|---|
| Protein embedding | ESM-family, ProtT5-family, or comparable protein language-model embeddings | Map amino-acid sequences into a learned protein landscape |
| Structure augmentation | AlphaFold/ESMFold-style predicted structures, confidence metrics, surface exposure, charge and hydrophobic patches | Add physical geometry when structure quality and rights allow |
| Condition augmentation | pH, temperature, concentration, buffer, storage, agitation, and stress metadata | Place the same protein under different physical contexts |
| Risk-direction learning | Linear probes, regularized models, gradient-boosted trees, metric learning, contrastive learning | Learn directions associated with aggregation, instability, low expression, poor SEC, or other endpoints |
| Danger-zone rules | Failure-anchor retrieval, clustering, calibrated thresholds, conformal or distance-based applicability checks | Convert embedding-space signals into reliable, reviewable risk regions |

The product insight is that embedding space alone is not the asset. The asset is a protein landscape annotated by experimentally grounded physical-risk directions.

### B. Lab-in-the-loop physical AI autopilot

Goal: turn the risk landscape into a closed-loop system that recommends what to test next, records frozen predictions, learns from wet-lab outcomes, and updates the reliability map.

```text
current risk landscape
    -> next experiment policy
    -> wet-lab / robot / partner assay
    -> frozen evidence registry
    -> calibration update
    -> next-round candidate selection
```

Key technical stack:

| Component | AI method | Technical role |
|---|---|---|
| Experiment selection | Active learning, uncertainty sampling, diversity-aware batch selection | Choose the most informative next candidate panel |
| Condition navigation | Bayesian experimental design or Bayesian optimization | Explore formulation and stress conditions efficiently |
| Evidence update | Calibration models, conformal prediction, uncertainty updates, error analysis | Convert wet-lab readouts into improved risk estimates |
| Workflow agents | Assay-planning agent, evidence-tracking agent, report-generation agent | Coordinate repeatable scientific operations without replacing human review |
| Automation connectors | ELN/LIMS integration, robot-ready assay manifests, API workflows | Connect prediction, experiment execution, and evidence capture |

The product insight is that AI agents should not be presented as autonomous scientists. They should be evidence-loop operators that keep predictions, experiments, and model updates aligned.

## Key technical questions

PRIB-KI should be judged by whether it can answer these questions better than simple heuristics:

1. Where is a candidate located in the protein foundation-model landscape?
2. Which physical failure anchors or danger zones is it near?
3. Does the learned risk direction improve over basic descriptors such as length, pI, charge, GRAVY, and hydrophobicity?
4. Is the candidate inside the validated applicability domain?
5. Which wet-lab readout would most reduce uncertainty?
6. Does a frozen prediction enrich real experimental failures in a prospective panel?
7. Can a partner process a candidate batch and trace every output to input data, model version, and assay context?

## Metrics that matter

The project should avoid vague AI claims and measure progress through explicit technical and product metrics.

| Layer | Metrics |
|---|---|
| Landscape quality | embedding coverage, cluster stability, sequence-family leakage checks, failure-anchor separation |
| Risk prediction | ROC-AUC/PR-AUC for binary endpoints, Spearman correlation for continuous endpoints, top-k enrichment for high-risk candidates |
| Baseline advantage | improvement over length, pI, charge, GRAVY, hydrophobicity, linear scores, and regularized baselines |
| Reliability | calibration error, conformal coverage, uncertainty ranking quality, out-of-domain detection rate |
| Prospective evidence | enrichment of predefined failures in high-risk groups, false-negative review, protocol deviation log |
| Product value | reduction in low-value wet-lab screening, time from candidate batch to ranked report, percentage of outputs with traceable evidence |

## PRIB-KI AI stack

The target architecture should be modular, so each layer can be tested independently.

```text
1. Data layer
   public sequences, partner candidates, assay endpoints, condition metadata

2. Reliability landscape layer
   sequence descriptors, protein language-model embeddings, predicted structures,
   geometric features, confidence metrics, perturbation descriptors, failure anchors

3. Physical reliability layer
   aggregation, stability, solubility, expression, viscosity, chemical liability,
   scale-up and handling-stress risk directions

4. Model layer
   interpretable baselines, calibrated ML models, uncertainty estimation,
   applicability-domain detection, endpoint-specific models, active learning

5. Decision layer
   candidate ranking, accept/review/reject grouping, assay planning,
   audit-ready reports and API output

6. Learning layer
   retrospective benchmarks, prospective wet-lab panels, partner feedback,
   model and threshold updates under version control
```

## Current implementation versus target system

| Layer | Current public prototype | Target system |
|---|---|---|
| Data | Public therapeutic-antibody sequences plus synthetic stress variants for demo behavior. | Endpoint-specific experimental datasets with assay conditions, protein-class labels, and partner panels. |
| Representation | Sequence-derived physicochemical descriptors and motif proxies. | Protein language-model landscape with structure/condition augmentation and failure anchors. |
| Physical logic | Interpretable heuristic links to aggregation, scale-up sensitivity, and stability axes. | Learned physical-risk directions and danger-zone rules calibrated to measured endpoints such as expression yield, SEC monomer percentage, DSF/nanoDSF, aggregation onset, solubility, viscosity, and self-interaction. |
| AI model | AI-weighted pseudo-target mapping for demonstrator behavior. | Endpoint-specific supervised models, baseline comparisons, uncertainty estimates, applicability-domain checks, active learning, and prospective calibration. |
| Product output | Ranked table, risk axes, wet-lab effort index, visualization, CSV export. | API, audit-ready reports, project-specific settings, assay recommendations, evidence registry, and model-version traceability. |

## Technical route

### Phase 0 - transparent demonstrator

Goal: prove the workflow is inspectable and reproducible.

Deliverables:

- deterministic demo dataset and fixed seeds;
- documented descriptor and scoring logic;
- clear claims boundary;
- unit tests for feature extraction and scoring;
- output schema for ranked candidate tables.

Evidence gate:

- a clean checkout reproduces the same ranking and the same configuration metadata.

### Phase 1 - benchmark-ready AI foundation

Goal: move from demo scoring to evidence-generating infrastructure.

Deliverables:

- benchmark data schema: sequence, protein class, assay endpoint, assay condition, readout, unit, source, usage rights;
- sequence-cluster, scaffold-aware, family-aware, or temporal split protocol;
- simple baselines: length, pI, charge, GRAVY, hydrophobicity, unweighted linear score, regularized linear model;
- reproducible benchmark scripts and confidence intervals;
- model-card style reporting for each endpoint.

Evidence gate:

- at least one measured endpoint is evaluated on an independent holdout set and compared with simple baselines.

### Phase 2 - protein physical reliability landscape

Goal: build the Stage A AI system: a protein foundation-model landscape annotated by physical-risk directions and danger-zone rules.

Candidate inputs:

- protein language-model embeddings for sequence context;
- predicted structures from suitable tools where usage rights and throughput allow;
- structure confidence metrics such as local confidence and pairwise error;
- surface exposure, charge patches, hydrophobic patches, interface features, and geometric descriptors;
- perturbation metadata such as pH, temperature, concentration, agitation, storage, and formulation context.

Modeling tasks:

- learn physical-risk directions for aggregation, stability, expression, SEC, solubility, or self-interaction endpoints;
- retrieve nearest stable and failed anchors in embedding space;
- separate reliable, fragile, and unknown regions through calibrated thresholds;
- report applicability domain based on sequence family, embedding distance, protein class, and missing context.

Design rule:

- new representations should be retained only when they improve endpoint performance under leakage-aware splits and remain interpretable enough for wet-lab review.

Evidence gate:

- the embedding-landscape model improves top-k failure enrichment or endpoint correlation over physicochemical baselines under family-aware splits.

### Phase 3 - calibrated physical-AI modules

Goal: convert PRIB-KI from relative screening into endpoint-specific risk modules.

Model families to compare:

- interpretable generalized linear or additive models;
- tree-based models for tabular descriptors;
- graph or geometric neural networks when structure data is reliable;
- multi-task models when endpoints share enough data;
- conformal prediction or calibrated uncertainty wrappers;
- applicability-domain models based on sequence similarity, descriptor distance, protein class, and missing context.

Evidence gate:

- the model improves over baselines on a held-out endpoint, reports uncertainty, flags out-of-domain candidates, and passes ablation and sensitivity analysis.

### Phase 4 - lab-in-the-loop physical AI

Goal: build the Stage B AI system: a closed evidence loop that uses frozen predictions and wet-lab readouts to update the risk landscape.

Deliverables:

- frozen code, model, thresholds, and configuration before wet-lab results are known;
- 30-50 candidate panel for the first prospective study;
- at least two orthogonal readouts, for example expression yield plus SEC monomer percentage or DSF/nanoDSF;
- prediction registry;
- failure analysis for both correct and incorrect predictions.
- active-learning policy for selecting a next candidate panel;
- assay recommendation logic based on uncertainty and missing evidence;
- model calibration update after readouts are received.

Evidence gate:

- high-risk candidates are enriched for predefined experimental problems compared with baselines, and all protocol deviations are recorded.

### Phase 5 - pilot product integration

Goal: turn the validated workflow into partner-usable software.

Deliverables:

- batch API and project-level configuration;
- audit-ready PDF/HTML reports;
- ELN/LIMS/export integration;
- robot-ready or CRO-ready assay manifests where partner infrastructure supports them;
- role-based access and data-handling rules;
- monitoring for model drift and endpoint-specific performance;
- partner-specific calibration while preserving a shared core evidence standard.

Evidence gate:

- a partner can process a candidate panel without manual code changes and trace every output to input data, model version, configuration, and assay context.

## Where PRIB-KI sits in the broader AI landscape

PRIB-KI should be described as downstream and complementary to major AI-biology systems:

| AI-biology capability | Examples of field direction | PRIB-KI relationship |
|---|---|---|
| Structure prediction | AlphaFold-family and ESMFold-style models predict structures and complexes. | Use predictions and confidence metrics as optional inputs, not as final developability evidence. |
| Generative protein design | Diffusion and inverse-folding methods can generate backbones and sequences. | Screen and prioritize generated candidates before costly experiments. |
| Protein language models | Sequence-scale foundation models learn representations useful for structure, function, or design tasks. | Use embeddings as candidate features after benchmark validation. |
| Experimental AI loops | Design-test-learn cycles link models with measured assay results. | Build the feedback loop around physical reliability endpoints. |

## Practical product modules

The platform should eventually expose the AI framework as modules that map to user decisions:

| Module | User question | Likely data inputs | Output |
|---|---|---|---|
| Reliability screen | Which candidates should enter wet-lab validation first? | Sequence, protein class, descriptors, optional structure | Ranked list and accept/review/reject grouping |
| Failure-mode explainer | Why is this candidate risky? | Feature contributions and physical modules | Human-readable risk drivers |
| Assay planner | Which experiment should we run next? | Uncertainty, missing evidence, risk profile | Suggested validation readouts |
| Applicability-domain guard | Can the model be trusted for this candidate? | Similarity, descriptor distance, class metadata | In-domain, review, or out-of-domain flag |
| Feedback calibrator | Did experiments confirm the ranking? | Frozen predictions and assay results | Model update package and validation report |

## Core technical principles

1. Keep physics-linked descriptors visible even when deeper AI models are added.
2. Treat sequence-only, structure-augmented, and endpoint-calibrated models as separate releases.
3. Never train or evaluate with random row splits when related protein families can leak across splits.
4. Compare every AI model against simple physicochemical baselines.
5. Report uncertainty and applicability domain before reporting a confident score.
6. Use synthetic variants only for software tests, never as experimental negative labels.
7. Let wet-lab results update the model, but freeze predictions before prospective validation.

## Reference landscape

The technical route is aligned with the direction of the field, where structure prediction, generative design, protein language models, and accessible folding workflows are becoming infrastructure for biological engineering:

- Abramson et al., "Accurate structure prediction of biomolecular interactions with AlphaFold 3," Nature, 2024. https://www.nature.com/articles/s41586-024-07487-w
- Watson et al., "De novo design of protein structure and function with RFdiffusion," Nature, 2023. https://www.nature.com/articles/s41586-023-06415-8
- Dauparas et al., "Robust deep learning-based protein sequence design using ProteinMPNN," Science, 2022. https://www.science.org/doi/10.1126/science.add2187
- Lin et al., "Evolutionary-scale prediction of atomic-level protein structure with a language model," Science, 2023. https://www.science.org/doi/10.1126/science.ade2574
- Mirdita et al., "ColabFold: making protein folding accessible to all," Nature Methods, 2022. https://www.nature.com/articles/s41592-022-01488-1
