# PRIB-KI

Failure-oriented protein developability screening.

Maintained by **PRIB-KI Lab**.

Try at:
https://prib-ki.streamlit.app/

Intro video (03:36):
https://youtu.be/2hRn8kXBsRI

## Status

PRIB-KI is a research prototype for early-stage protein candidate prioritization. The current implementation focuses on sequence-derived descriptors, heuristic failure-mode features, relative risk aggregation, and interactive batch review.

The current release is not an experimentally validated predictor of clinical or manufacturing failure. Scores are relative screening indices and should be used to support, not replace, wet-lab developability testing.

## Current capabilities

- protein sequence validation;
- physicochemical descriptor calculation;
- motif-based liability features;
- failure-oriented risk axes;
- configurable score aggregation;
- batch ranking and filtering;
- Streamlit-based review interface;
- deterministic demonstration dataset generation.

## Repository layout

```text
.
├── app.py                  # Streamlit application
├── risk_engine.py          # feature and scoring implementation
├── download_data.py        # source data retrieval
├── build_demo_set.py       # deterministic demo dataset builder
├── data/                   # demonstration data
└── docs/
    ├── methodology.md      # scientific and computational assumptions
    ├── validation.md       # benchmark and prospective validation plan
    ├── roadmap.md          # engineering milestones
    └── data.md             # data provenance and usage limits
```

## Method summary

```text
sequence input
    -> input validation
    -> descriptor calculation
    -> failure-mode features
    -> risk aggregation
    -> candidate ranking
    -> review / retain / deprioritize
```

The project follows a failure-oriented approach: a plausible structure is not assumed to imply robustness under expression, purification, formulation, storage, concentration, agitation, or other non-ideal conditions.

See [`docs/methodology.md`](docs/methodology.md) for the current assumptions and implementation boundaries.

## Installation

Python 3.12 is recommended.

```bash
python -m venv .venv
```

Activate the environment and install the dependencies:

```bash
pip install -r requirements.txt
```

## Demo data

Download the public source file and build the deterministic demo set:

```bash
python download_data.py
python build_demo_set.py
```

The generated demo contains public therapeutic-antibody sequences and synthetic stress variants. Synthetic variants are software test cases, not experimental negative labels. See [`docs/data.md`](docs/data.md).

## Run

```bash
streamlit run app.py
```

## Validation status

The repository currently demonstrates workflow execution and controlled descriptor stratification.

Before quantitative performance claims are made, the system requires:

- endpoint-specific experimental benchmark data;
- sequence- or family-aware train/test separation;
- comparison with simple physicochemical baselines;
- uncertainty and applicability-domain analysis;
- blinded prospective validation.

The planned evaluation procedure is documented in [`docs/validation.md`](docs/validation.md).

## Intended use

Appropriate uses:

- research prototyping;
- early candidate triage;
- descriptor and hypothesis exploration;
- workflow demonstrations;
- planning of experimental validation.

Out-of-scope uses:

- clinical decisions;
- regulatory submissions;
- manufacturing release decisions;
- absolute probability of development failure;
- replacement of experimental assays.

## Development roadmap

See [`docs/roadmap.md`](docs/roadmap.md).

## References

- Chen, L. T. et al. Target sequence-conditioned design of peptide binders using masked language modeling. *Nature Biotechnology* (2025).
- Sappington, I. et al. Improved protein binder design using beta-pairing targeted RFdiffusion. *Nature Communications* (2026).
- Didi, K. et al. Scaling Atomistic Protein Binder Design with Generative Pretraining and Test-Time Compute. arXiv:2603.27950 (2026).

## Contact

PRIB-KI Lab  
Dr. Haipei Liu  
haipei.thu@gmail.com
