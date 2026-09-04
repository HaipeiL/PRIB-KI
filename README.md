<p align="center">
  <img src="logo.png" alt="PRIB-KI logo" width="220">
</p>

# PRIB-KI

**Physical-AI reliability screening for AI-designed proteins.**

PRIB-KI is an early-stage research and commercialization project building a **failure-first reliability layer** between computational protein design and expensive wet-lab development. The central question is not only whether a protein can be generated or folded, but whether it is **stable and reliable enough to deserve experimental investment**.

[Live demonstrator](https://prib-ki.streamlit.app/) | [September 2026 R&D progress](docs/research_progress_2026.md) | [AI framework](docs/ai_framework.md) | [Methodology](docs/methodology.md) | [Validation](docs/validation.md) | [Long-term vision](docs/long_term_vision.md) | [Repository review guide](docs/repository_review_guide.md)

---

## September 2026 R&D snapshot

PRIB-KI has progressed beyond the original sequence-only demonstrator. Recent work combines transparent physicochemical baselines with protein foundation models, perturbation-response learning, historical wet-lab workflow data, and multi-GPU execution.

A completed parameter-efficient adaptation of **ESM-2 3B** on experimental protein-stability data achieved the following on the implemented held-out validation split with **30,000 evaluation examples**:

| Evaluation target | Metric | Result |
|---|---|---:|
| Absolute protein stability | Spearman | **0.8875** |
| Absolute protein stability | Pearson | **0.8597** |
| Perturbation-induced stability change | Spearman | **0.8276** |
| Perturbation-induced stability change | Pearson | **0.8237** |
| Direction of stability change | Sign accuracy | **82.53%** |

The same R&D program has also validated distributed ESM-2 3B perturbation scoring at practical throughput: a retained run scored **30,400 sequence mutants across 100 proteins in ~4.4 minutes on two H100 GPUs**.

These are **held-out research results under the implemented dataset split**, not yet a claim of universal generalization to unseen protein families or industrial conditions. Family-aware/external holdouts and prospective wet-lab validation remain explicit next gates. See [`docs/research_progress_2026.md`](docs/research_progress_2026.md) for the full public technical summary and claim boundaries.

## Why PRIB-KI exists

Modern protein AI can generate, rank, and structurally evaluate more candidates than most experimental teams can afford to test. This moves an important bottleneck downstream:

- Which candidates should enter limited wet-lab validation capacity?
- Which designs carry hidden stability or developability risk despite plausible sequence or structure?
- Which perturbations are likely to increase or decrease reliability?
- Which experimental measurement would most reduce uncertainty before further investment?

PRIB-KI addresses this gap as a reliability and decision-support layer:

```text
AI protein design / candidate libraries
        -> PRIB-KI reliability assessment
        -> prioritized wet-lab validation
        -> experimental evidence and calibration
        -> better next-round selection
```

The product thesis is intentionally narrow: **reduce low-value experimental screening by identifying fragile or uncertain candidates earlier and making the reason for prioritization inspectable.**

## What makes the approach different

PRIB-KI is not positioned as another protein generator or structure predictor. Its scientific framing is **failure-first** rather than success-first.

```text
perturbation
    -> molecular response
    -> failure mechanism
    -> experimental readout
    -> computational representation
    -> decision signal
```

The technical program combines four layers:

1. **Transparent physical baselines** — sequence descriptors, charge, hydrophobicity, motif and liability signals remain visible and auditable.
2. **Protein foundation-model representations** — learned sequence representations are adapted to stability and perturbation-response tasks.
3. **Experimental evidence infrastructure** — historical and future partner readouts are organized around explicit stages, endpoints, conditions, and failure definitions.
4. **Decision and feedback layer** — candidate ranking, applicability-domain checks, uncertainty, and ultimately frozen prospective predictions followed by wet-lab feedback.

The goal is not a universal black-box "protein quality score." The goal is a traceable **protein physical reliability landscape** that can support go/no-go and next-experiment decisions.

## Current technical status

| Layer | Public demonstrator | Active R&D status | Next validation gate |
|---|---|---|---|
| Input | Sequence CSV (`id`, `sequence`) | Experimental stability records and perturbation pairs | Condition- and protein-class-aware manifests |
| Representation | Physicochemical descriptors and motif proxies | ESM-2 3B protein-language-model representations | Leakage-aware family/scaffold/external holdouts |
| Stability learning | Heuristic demo scoring | Completed supervised stability and perturbation-response adaptation | Baseline comparison and external/family-aware validation |
| Failure evidence | TargetTrack historical workflow integration | Large-scale experimental-history processing | Endpoint-specific industrial/wet-lab datasets |
| Output | Risk axes, ranking, effort proxy, CSV export | Learned stability / perturbation signals and research evaluation | Calibrated uncertainty and applicability domain |
| Product loop | Interactive demonstrator | Reproducible model/HPC workflows | Frozen prospective predictions + wet-lab feedback |

The public Streamlit application is intentionally simpler than the active research stack. It demonstrates the user workflow and interpretable risk logic without releasing unpublished model assets, training data, detailed split assignments, or commercialization-sensitive implementation choices.

## Historical wet-lab workflow integration

PRIB-KI has completed a reproducible integration of the Protein Structure Initiative's final 2017 TargetTrack archive.

Verified local snapshot:

| Unit | Records |
|---|---:|
| Protein targets | **335,771** |
| Trial / construct units | **961,548** |
| Historical status events | **3,783,070** |
| Unmapped source statuses | **0** |

The pipeline reconstructs experimental progression such as:

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

It provides an auditable schema for trial-level progression, repeated attempts, stop reasons, and target-level aggregation. TargetTrack is used as historical workflow/evidence infrastructure, **not as a direct predictive-validation dataset for the public demo**. See [`docs/targettrack_wetlab_failure_funnel.md`](docs/targettrack_wetlab_failure_funnel.md).

## Public demonstrator

The current Streamlit demonstrator supports:

- CSV input with `id` and `sequence` columns;
- sequence validation and interpretable physicochemical descriptors;
- motif-based liability features;
- nine sub-risk signals;
- three decision-oriented axes: `Aggregation`, `ScaleUpSensitivity`, and `Stability`;
- configurable risk weighting;
- a wet-lab effort proxy;
- accept / review / reject grouping;
- batch ranking, visualization, filtering, and CSV export;
- single-candidate explainability.

The public scoring engine is implemented in [`risk_engine.py`](risk_engine.py), and the application in [`app.py`](app.py).

### Demo data

The bundled demonstration set contains:

- 60 public therapeutic-antibody sequences sampled from Thera-SAbDab;
- 40 deterministic synthetic stress variants;
- fixed seed `42` for reproducibility.

Synthetic variants are software test cases, **not experimentally failed proteins** and not negative labels for predictive-performance reporting. See [`docs/data.md`](docs/data.md).

## Research execution capability

PRIB-KI maintains reproducible multi-GPU workflows for protein foundation-model inference, perturbation scoring, and parameter-efficient adaptation. Recent retained workflows have used institutional GPU/HPC infrastructure including **WestAI / RWTH Aachen** and **hessian.AI Cluster 42**.

The public repository does not expose infrastructure identifiers, internal environment runbooks, private model caches, checkpoints, or detailed job records. Those are maintained separately for reproducibility and operational security.

## AI-to-lab execution experience

The founding team has also documented prior AI-to-lab execution in the GEM x Adaptyv RBX1 Binder Design Competition. **7 of 21 submitted candidates advanced to wet-lab testing**, corresponding to approximately the top **2.5% of more than 12,000 global submissions**. No confirmed binder was ultimately obtained; the case is therefore used as execution evidence and as a practical motivation for stronger post-generation reliability screening rather than as PRIB-KI model validation.

See [`docs/protein_design_execution_case.md`](docs/protein_design_execution_case.md).

## Founding team

PRIB-KI combines protein biophysics, physics-informed computation, and commercialization rather than treating protein reliability as a generic software problem.

- **Dr. Haipei Liu — Founder & scientific lead.** Background in single-molecule biophysics, protein mechanics, AFM-based force spectroscopy, quantitative modeling, and AI-assisted scientific data analysis; research experience at the University of Hong Kong, University of Basel, and ETH Zurich/DBSSE.
- **Dr. Xianglong Peng — Co-founder & scientific lead.** Background in solid mechanics, nonlinear/energy-based modeling, computational mechanics, and deep-learning approaches to complex physical systems; currently affiliated with TU Darmstadt.
- **Bingyi Zhao — Co-founder, business development & market interface.** Background in industrial engineering and management, technology-oriented consulting, international project coordination, and commercialization strategy.

The team structure reflects the core PRIB-KI thesis: protein reliability sits at the intersection of **biophysical failure mechanisms, scalable AI, and real experimental decision-making**.

## Validation standard

PRIB-KI distinguishes between an encouraging internal research result and a validated product claim.

The current foundation-model result is a strong held-out benchmark under the implemented split. Before broader predictive or industrial claims, the validation program requires:

- sequence-cluster, family-aware, scaffold-aware, temporal, or external holdout testing where appropriate;
- comparison against simple physicochemical and regularized machine-learning baselines under identical splits;
- uncertainty and applicability-domain analysis;
- endpoint-specific reporting rather than an opaque composite score;
- frozen prospective predictions before wet-lab results are known;
- explicit analysis of false negatives and protocol deviations.

A practical prospective study is a **30–50 candidate panel** with at least two orthogonal experimental readouts.

See [`docs/validation.md`](docs/validation.md).

## What PRIB-KI is not

PRIB-KI is not:

- a protein generator;
- a replacement for AlphaFold or other structure-prediction systems;
- a replacement for developability assays;
- a clinical or regulatory decision tool;
- a manufacturing release system;
- a claim that current internal validation already proves performance across all protein families or industrial conditions.

The intended role is **pre-wet-lab prioritization and evidence-guided decision support**.

## Repository map

```text
.
├── app.py                         # Streamlit public demonstrator
├── risk_engine.py                 # transparent descriptor/scoring baseline
├── build_demo_set.py              # deterministic demo builder
├── data/                          # public/demo data only
├── prib_targettrack/              # reproducible TargetTrack processing
├── docs/
│   ├── research_progress_2026.md  # current quantitative R&D snapshot
│   ├── ai_framework.md            # AI architecture and technical route
│   ├── methodology.md             # scientific assumptions and boundaries
│   ├── validation.md              # evidence standards and next gates
│   ├── targettrack_wetlab_failure_funnel.md
│   ├── protein_design_execution_case.md
│   ├── long_term_vision.md
│   ├── roadmap.md
│   ├── data.md
│   └── repository_review_guide.md
└── requirements.txt
```

## 5-minute review path

1. Read the **September 2026 R&D snapshot** above.
2. Open the [live demonstrator](https://prib-ki.streamlit.app/) and run the built-in dataset.
3. Read [`docs/research_progress_2026.md`](docs/research_progress_2026.md) for the completed foundation-model results and disclosure boundaries.
4. Read [`docs/ai_framework.md`](docs/ai_framework.md) to understand the protein reliability landscape and lab-in-the-loop route.
5. Read [`docs/validation.md`](docs/validation.md) to see exactly what PRIB-KI still requires before stronger generalization or product claims.

## Local installation

Python 3.12 is recommended for the public demonstrator.

```bash
python -m venv .venv
pip install -r requirements.txt
streamlit run app.py
```

Optional demo-data rebuild:

```bash
python download_data.py
python build_demo_set.py
```

## Public disclosure boundary

PRIB-KI is simultaneously an active research program and an early commercialization project. This repository therefore publishes enough information to make the scientific direction, present evidence, data provenance, and product logic inspectable, while intentionally excluding:

- unpublished model checkpoints;
- private or restricted research datasets;
- exact internal training recipes and search history;
- detailed split assignments intended for publication-quality analysis;
- partner data and partner-specific calibration;
- internal HPC credentials, paths, job records, and operating runbooks;
- commercialization-sensitive implementation details.

## License and reuse

No open-source license has been declared yet. Until a license file is added, all rights are reserved by default. Please contact PRIB-KI Lab before reuse, redistribution, or commercial evaluation outside the intended demonstration context.

## Contact

**PRIB-KI Lab**  
Dr. Haipei Liu  
haipei.thu@gmail.com
