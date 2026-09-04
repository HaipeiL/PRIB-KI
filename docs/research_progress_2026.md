# PRIB-KI Research Progress — September 2026

Maintained by **PRIB-KI Lab**.

This page summarizes the most important technical progress that can be stated publicly without releasing unpublished model details, research data, partner information, or implementation choices that remain part of ongoing publication and commercialization work.

## Current status in one sentence

PRIB-KI has progressed from a transparent sequence-level demonstrator to an active protein-foundation-model R&D program with experimentally grounded stability learning, high-throughput perturbation scoring, and reproducible historical wet-lab data integration.

## Completed R&D highlights

### 1. Protein language-model adaptation for stability and perturbation response

In August 2026, PRIB-KI completed a parameter-efficient adaptation of an **ESM-2 3B** protein language model using experimental protein-stability labels.

On the implemented held-out validation split, evaluated on **30,000 examples**, the model achieved:

| Evaluation target | Metric | Result |
|---|---|---:|
| Absolute protein stability | Spearman correlation | **0.8875** |
| Absolute protein stability | Pearson correlation | **0.8597** |
| Perturbation-induced stability change | Spearman correlation | **0.8276** |
| Perturbation-induced stability change | Pearson correlation | **0.8237** |
| Direction of stability change | Sign accuracy | **82.53%** |
| Learned embedding displacement vs. magnitude of stability change | Spearman correlation | **0.5136** |

These results support three practical conclusions under the current validation protocol:

- the learned representation carries strong signal for ranking absolute stability;
- it captures the direction and relative magnitude of sequence-perturbation effects;
- stability direction is currently learned more strongly than exact geometric distance in embedding space.

**Important boundary:** these numbers are held-out validation results under the implemented dataset split. They are not yet a claim of generalization to unseen protein families, industrial conditions, or prospective partner datasets. Leakage-aware family/scaffold holdouts and prospective validation remain explicit next gates.

### 2. High-throughput perturbation scoring

PRIB-KI has also validated a distributed protein-language-model perturbation workflow for rapid sequence-level scoring.

A retained ESM-2 3B execution processed:

- **100 proteins**;
- **30,400 sequence mutants**;
- on **2 NVIDIA H100 GPUs**;
- in approximately **4.4 minutes** end-to-end for the retained run.

This result is primarily an execution and scalability demonstration rather than a biological-performance benchmark. It shows that foundation-model perturbation scoring can be integrated into a practical candidate-screening workflow without requiring very large GPU fleets.

### 3. Historical wet-lab execution data integration

The public repository already includes the reproducible TargetTrack historical workflow. The verified local snapshot contains:

| Unit | Records |
|---|---:|
| Protein targets | **335,771** |
| Trial / construct units | **961,548** |
| Historical status events | **3,783,070** |
| Unmapped source statuses | **0** |

This work provides PRIB-KI with an auditable schema for experimental-stage progression, stop reasons, repeated attempts, and trial-versus-target aggregation. It is used as workflow and evidence-infrastructure research, not as a direct performance label set for the public risk demonstrator.

## How the technical program now fits together

```text
protein sequence / candidate library
        |
        +--> transparent descriptor baseline
        |
        +--> protein language-model representation
                    |
                    +--> absolute stability signal
                    +--> perturbation-response signal
                    +--> learned reliability directions
        |
        +--> historical / partner experimental evidence
                    |
                    +--> failure-stage taxonomy
                    +--> endpoint-specific calibration
                    +--> applicability-domain checks
        |
        v
candidate prioritization before expensive wet-lab work
```

The public Streamlit application remains intentionally simpler than the active research stack. Its purpose is to make the decision workflow, risk dimensions, assumptions, and interface inspectable without publishing research assets that are still under development.

## Public repository vs. active research stack

| Area | Public repository | Active R&D |
|---|---|---|
| Sequence descriptors | Included | Included as transparent baselines |
| Demo risk workflow | Included | Used as product/interface reference |
| Protein language-model work | Architecture and results summarized | Active training and validation |
| Model checkpoints | Not released | Retained internally |
| Experimental training data | Not redistributed here | Used under applicable source/data constraints |
| Exact split assignments and preprocessing | Not released in full | Retained for reproducibility and publication work |
| Detailed hyperparameter search | Not released | Internal R&D record |
| Partner data | Never public by default | Isolated under partner-specific governance |
| Production/HPC operating details | Not public | Internal reproducibility records |

This separation is deliberate. PRIB-KI is being developed simultaneously as a research program and a commercialization project; the repository is intended to demonstrate technical credibility and scientific discipline without functioning as a full release of unpublished IP.

## What the results do — and do not — establish

The completed work establishes that PRIB-KI can:

- build and operate reproducible protein-foundation-model workflows at multi-GPU scale;
- learn strong held-out stability and perturbation-response signals from experimental labels;
- connect model development to an explicit failure-oriented product thesis;
- maintain transparent baselines and public claim boundaries while deeper models remain under active research;
- structure large historical experimental records for future evidence-loop integration.

The completed work does **not** yet establish:

- universal protein-family generalization;
- prospective industrial developability performance;
- calibrated probabilities of aggregation, formulation failure, or manufacturing failure;
- clinical or regulatory validity;
- replacement of wet-lab assays.

## Next evidence gates

The next R&D gates are deliberately stricter than the current internal holdout result:

1. sequence-cluster / family-aware or external holdout evaluation;
2. comparison against physicochemical and simpler machine-learning baselines under identical splits;
3. endpoint-specific uncertainty and applicability-domain analysis;
4. frozen prospective predictions on a new experimental panel;
5. integration of measured wet-lab feedback into the reliability landscape.

The target is not a black-box "protein quality score." The target is a traceable reliability layer that helps experimental teams decide **which candidates are worth testing next, why, and with what level of confidence**.

## Disclosure policy

PRIB-KI publishes enough information to make the scientific direction, current evidence, and product logic inspectable. Exact unpublished training recipes, private datasets, partner information, model checkpoints, infrastructure identifiers, and commercialization-sensitive implementation details are intentionally excluded until they are appropriate to release.
