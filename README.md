<p align="center">
  <img src="logo.png" alt="PRIB-KI logo" width="220">
</p>

# PRIB-KI

**Physical-AI reliability screening for AI-designed proteins.**

PRIB-KI is a research prototype for early-stage protein candidate triage. It helps protein design, biotech, pharma, and enzyme-engineering teams decide which candidates deserve wet-lab budget before fragile designs consume expression, formulation, stability, or scale-up experiments.

[Live demonstrator](https://prib-ki.streamlit.app/) | [Intro video](https://youtu.be/2hRn8kXBsRI) | [Long-term vision](docs/long_term_vision.md) | [Execution case](docs/protein_design_execution_case.md) | [AI framework](docs/ai_framework.md) | [TargetTrack foundation](docs/targettrack_wetlab_failure_funnel.md) | [Methodology](docs/methodology.md) | [Validation plan](docs/validation.md) | [Roadmap](docs/roadmap.md) | [Funder and investor brief](docs/funder_investor_brief.md) | [Repository review guide](docs/repository_review_guide.md)

---

## 5-minute review path

1. Open the [live demonstrator](https://prib-ki.streamlit.app/) and run the built-in demo dataset.
2. Review the implemented scoring logic in [`risk_engine.py`](risk_engine.py).
3. Read the concise [`funder_investor_brief.md`](docs/funder_investor_brief.md) and the validation boundaries in [`validation.md`](docs/validation.md).

## Why PRIB-KI exists

Generative protein design and structure prediction can produce large numbers of plausible protein candidates. The downstream bottleneck is increasingly not generation, but selection:

- Which candidates should enter limited wet-lab validation capacity?
- Which designs carry hidden developability risk despite plausible sequence or structure?
- Which candidates are likely to fail under formulation, storage, concentration, handling, production, or scale-up stress?

PRIB-KI addresses this gap as a **failure-first reliability layer** between in-silico design and experimental validation.

```text
AI protein design / candidate libraries
        -> PRIB-KI reliability screening
        -> prioritized wet-lab validation
        -> experimental calibration and feedback
```

The current public release focuses on sequence-derived descriptors and interpretable risk signals. It is designed to demonstrate workflow feasibility, transparent assumptions, and a path toward experimental validation.

## Long-term vision

The company-level vision is documented in [`docs/long_term_vision.md`](docs/long_term_vision.md). In short, PRIB-KI is built on the judgement that AI protein engineering is moving from a generation bottleneck to a reliability bottleneck.

The long-term route has two parts: build a **protein physical reliability landscape**, then use **lab-in-the-loop physical AI** to navigate and update that landscape through wet-lab evidence.

The team has documented prior AI-to-lab execution in the GEM x Adaptyv RBX1 Binder Design Competition, where 7 of 21 submitted candidates advanced to wet-lab testing, corresponding to approximately the top 2.5% of more than 12,000 global submissions. No confirmed binder was ultimately obtained, which reinforces PRIB-KI's focus on post-generation reliability screening. See [`docs/protein_design_execution_case.md`](docs/protein_design_execution_case.md).

## AI framework in one view

PRIB-KI should be read as an **AI for life science** project with a narrow, testable role in the drug x AI workflow:

```text
generative protein AI / discovery libraries
        -> physical-AI reliability screen
        -> wet-lab assay planning
        -> experimental feedback and calibration
```

The platform is designed around three connected AI ideas:

- **Drug x AI:** PRIB-KI supports therapeutic-protein, antibody, enzyme, and engineered-binder teams at the candidate-selection step before expensive developability experiments.
- **Physical AI:** PRIB-KI combines learned representations with biophysical failure logic, so aggregation, instability, self-interaction, viscosity, chemical liability, formulation stress, and scale-up sensitivity remain interpretable.
- **AI for life science:** PRIB-KI is intended to become a design-test-learn layer that links sequence, structure, foundation-model features, assay endpoints, and partner feedback.

See [`docs/ai_framework.md`](docs/ai_framework.md) for the AI architecture and technical route.

## What the prototype does

The Streamlit demonstrator supports:

- CSV input with `id` and `sequence` columns;
- protein sequence validation;
- physicochemical descriptor calculation;
- motif-based liability features;
- nine interpretable sub-risk signals;
- three decision-oriented risk axes:
  - `Aggregation`;
  - `ScaleUpSensitivity`;
  - `Stability`;
- configurable weighted aggregation;
- wet-lab effort index as a planning proxy;
- accept / review / reject grouping;
- batch ranking, filtering, visualization, and CSV export;
- single-candidate explainability with radar plots and sequence export.

## What PRIB-KI is not

The current repository is intentionally scoped and should not be read as a clinical, regulatory, or manufacturing decision tool.

PRIB-KI is not:

- a protein generator;
- a structure predictor;
- a replacement for developability assays;
- a calibrated probability model for industrial or clinical failure;
- a regulatory submission or manufacturing release system;
- a claim that synthetic demo variants are experimentally failed proteins.

Scores in this release are **relative screening indices** for research prototyping and workflow discussion. Experimental benchmark data and prospective wet-lab validation are required before quantitative performance claims are made.

## Current implementation

| Layer | Implemented now | Planned next |
|---|---|---|
| Input | Sequence CSV with `id` and `sequence` | Condition metadata, protein-class metadata, partner-specific manifests |
| Representation | Length, MW, pI, charge, GRAVY, hydrophobicity, motif proxies | Protein language-model embeddings, structure-derived descriptors, perturbation metadata, and applicability-domain checks |
| Risk logic | Interpretable heuristic features and AI-weighted pseudo-target mapping for demo use | Endpoint-specific calibration against experimental data with uncertainty estimates |
| Output | Risk axes, overall ranking, effort proxy, accept/review/reject exports | Audit-ready reports, API output, assay suggestions, model-version traceability |
| Validation | Workflow execution and deterministic demo dataset | Retrospective benchmarks, baseline comparisons, and blinded prospective wet-lab panel |

The risk engine is implemented in [`risk_engine.py`](risk_engine.py). The web interface is implemented in [`app.py`](app.py).

## Repository map

```text
.
├── app.py                     # Streamlit application
├── risk_engine.py             # Feature extraction and scoring implementation
├── download_data.py           # Public source-data retrieval helper
├── build_demo_set.py          # Deterministic demo dataset builder
├── data/
│   ├── TheraSAbDab.csv        # Public source data used for demo construction
│   └── demo_100.csv           # Demonstration dataset
├── docs/
│   ├── methodology.md         # Scientific and computational assumptions
│   ├── long_term_vision.md    # Company-level vision and two-route strategy
│   ├── protein_design_execution_case.md
│   │                           # Prior AI-to-lab execution case
│   ├── ai_framework.md        # AI positioning, architecture, and technical route
│   ├── validation.md          # Benchmark and prospective validation plan
│   ├── roadmap.md             # Technical de-risking milestones
│   ├── data.md                # Data provenance and usage limits
│   ├── funder_investor_brief.md
│   └── repository_review_guide.md
└── requirements.txt
```

## Scientific framing

PRIB-KI follows a failure-oriented evidence chain:

```text
perturbation
    -> molecular response
    -> failure mechanism
    -> experimental readout
    -> computational proxy
    -> decision signal
```

The current sequence-only implementation approximates early risk signals such as hydrophobic patch tendency, near-neutral charge behavior, chemical-liability motifs, and simple processability proxies. These signals are useful for demonstrator-level triage and hypothesis generation, but they are not direct substitutes for endpoint-specific experimental assays.

The next scientific step is to connect risk modules to measured endpoints such as expression yield, SEC monomer percentage, nanoDSF/DSF stability, aggregation onset, solubility, self-interaction, or viscosity.

## Completed data integration: TargetTrack historical wet-lab funnel

PRIB-KI has completed its first reproducible historical wet-lab data
integration using the Protein Structure Initiative's final 2017
[TargetTrack archive](https://doi.org/10.5281/zenodo.821654). The official
archive is downloaded locally, checked against its published MD5, inspected,
and streamed into analysis-ready experimental-history records. Raw archives
and row-level derived data remain excluded from Git.

### Dataset scope

| Unit | Local verified snapshot |
|---|---:|
| protein targets | 335,771 |
| trial / construct units | 961,548 |
| historical status events | 3,783,070 |
| unmapped source statuses | 0 |

### What the funnel records

```text
selected
  -> cloned
  -> expressed
  -> soluble
  -> purified
  -> crystallized
  -> diffraction-quality crystals
  -> structure determined
  -> deposited
```

For every trial/construct and target, the pipeline can identify the highest
observed stage, later successful retries, inferred prerequisite stages, and
whether a record ended in an explicit technical failure, nontechnical stop,
unknown stop, or censoring. This provides a historical reference for the
sequence of experimental work, not merely a count of finished structures.

### What this contributes to PRIB

This integration establishes the data workflow that PRIB needs when connecting
AI candidate prioritisation to wet-lab execution:

```text
candidate selected
  -> experimental stage recorded
  -> outcome and stop reason classified
  -> historical or partner evidence reviewed
  -> next experiment and next design round prioritised
```

It supplies a reusable schema, provenance record, stage map, stop taxonomy,
trial-versus-target aggregation, and auditable funnel reporting. Future
partner datasets can use the same structure with richer fields such as assay
conditions, expression yield, SEC, solubility, stability and explicit failure
reasons.

### Important boundary

The TargetTrack snapshot has 173,989 terminal trial records labelled `work stopped`,
but no explicit technical-failure category in the parsed status
history. They are therefore retained as unknown outcomes; the code does not
turn them into molecular failures or report a resolved-failure probability.
This is a historical structural-genomics workflow reference, not a validated
industrial developability model and it does not alter the Streamlit
demonstrator or its rankings.

See [the TargetTrack workflow](docs/targettrack_wetlab_failure_funnel.md) and
[the attribution notice](THIRD_PARTY_DATA.md).

## Demo data

The bundled demonstration set contains:

- 60 public therapeutic-antibody sequences sampled from Thera-SAbDab;
- 40 synthetic stress variants generated with deterministic sequence perturbations;
- fixed seed `42` for reproducibility;
- neutral candidate IDs after shuffling.

Synthetic variants are software test cases. They are not negative labels and must not be used to report ROC-AUC, accuracy, or predictive performance.

See [`docs/data.md`](docs/data.md) for data provenance and usage limits.

## Installation

Python 3.12 is recommended.

```bash
python -m venv .venv
```

Activate the environment, then install dependencies:

```bash
pip install -r requirements.txt
```

Run the app:

```bash
streamlit run app.py
```

Optional: rebuild the deterministic demo set after downloading the source file:

```bash
python download_data.py
python build_demo_set.py
```

## Validation standard

Before PRIB-KI is presented as a predictive model, the validation plan requires:

- deterministic software behavior for fixed inputs and configuration;
- explicit descriptor tests and missing-value handling;
- endpoint-specific experimental benchmark datasets;
- sequence-cluster, family-aware, scaffold-aware, or external holdout splits;
- comparison against simple physicochemical baselines;
- uncertainty and applicability-domain analysis;
- blinded prospective validation with frozen predictions.

A practical first prospective study is a 30-50 candidate panel with at least two orthogonal readouts.

See [`docs/validation.md`](docs/validation.md) for the evaluation protocol.

## Development roadmap

The roadmap is organized around technical and translational risk reduction, with explicit metrics for whether the AI route is working:

- `v0.2` reproducible demo and basic tests;
- `v0.3` endpoint-specific benchmark build;
- `v0.4` protein embedding landscape and calibrated physical-risk directions;
- `v0.5` prospective validation with frozen predictions;
- `v0.6` lab-in-the-loop pilot integration with API, evidence registry, and assay recommendations.

Target evidence metrics include reproducible rankings, top-risk enrichment for experimental failures, improvement over physicochemical baselines, calibration quality, out-of-domain detection, and reduction of low-value wet-lab screening effort in partner panels.

See [`docs/roadmap.md`](docs/roadmap.md).

## For funders, partners, and investors

This repository is meant to make the project inspectable:

- the current prototype can be run and reviewed;
- scientific assumptions are documented;
- data limitations are explicit;
- validation requirements are defined before performance claims;
- the commercial wedge is narrow and testable: reduce low-value wet-lab preselection effort by prioritizing physically more reliable candidates earlier.

For a concise due-diligence view, see [`docs/funder_investor_brief.md`](docs/funder_investor_brief.md).

## Responsible scope

The public prototype evaluates candidate sequences for reliability-related risk signals. It does not generate new biological designs, optimize activity, select clinical candidates, or recommend release decisions. Any future partner workflow should include data governance, assay documentation, model-version tracking, and human scientific review.

## License and reuse

No open-source license has been declared yet. Until a license file is added, all rights are reserved by default. Please contact PRIB-KI Lab before reuse, redistribution, or commercial evaluation outside the intended demonstration context.

## Contact

PRIB-KI Lab  
Dr. Haipei Liu  
haipei.thu@gmail.com
