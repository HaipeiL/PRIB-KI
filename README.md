<p align="center">
  <img src="logo.png" alt="PRIB-KI logo" width="220">
</p>

# PRIB-KI

**Physical reliability screening for AI-designed proteins.**

PRIB-KI is a research prototype for early-stage protein candidate triage. It helps protein design, biotech, pharma, and enzyme-engineering teams decide which candidates deserve wet-lab budget before fragile designs consume expression, formulation, stability, or scale-up experiments.

[Live demonstrator](https://prib-ki.streamlit.app/) | [Intro video](https://youtu.be/2hRn8kXBsRI) | [Methodology](docs/methodology.md) | [Validation plan](docs/validation.md) | [Roadmap](docs/roadmap.md) | [Funder and investor brief](docs/funder_investor_brief.md) | [Repository review guide](docs/repository_review_guide.md)

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
| Representation | Length, MW, pI, charge, GRAVY, hydrophobicity, motif proxies | Structure-derived descriptors and applicability-domain checks |
| Risk logic | Interpretable heuristic features and AI-weighted pseudo-target mapping for demo use | Endpoint-specific calibration against experimental data |
| Output | Risk axes, overall ranking, effort proxy, accept/review/reject exports | Audit-ready reports, API output, uncertainty estimates |
| Validation | Workflow execution and deterministic demo dataset | Retrospective benchmarks and blinded prospective wet-lab panel |

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

The roadmap is organized around technical and translational risk reduction:

- `v0.2` reproducible demo and basic tests;
- `v0.3` endpoint-specific benchmark build;
- `v0.4` calibrated research prototype;
- `v0.5` prospective validation;
- `v0.6` pilot integration with API and audit-ready reports.

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
