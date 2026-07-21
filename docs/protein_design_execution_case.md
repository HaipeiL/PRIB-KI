# Protein Design Execution Case

Maintained by **PRIB-KI Lab**.

This case note documents prior execution experience that shaped PRIB-KI's product direction. It should be read as team execution evidence, not as validation evidence for the current PRIB-KI scoring model.

## Case summary

The team participated in the **GEM x Adaptyv RBX1 Binder Design Competition**, connected to the **ICLR 2026 GEM Workshop**. The task was to design de novo protein binders against **RBX1**, a component of the Cullin-RING E3 ubiquitin ligase system.

The competition received more than **12,000 global design submissions**. From the PRIB-associated submitted set, **7 of 21 candidates advanced to the experimental stage**, corresponding to approximately the **top 2.5% of global submissions** by external selection performance.

The selected candidates were moved into the requested wet-lab expression and binding-evaluation workflow. The campaign did **not** ultimately produce a confirmed functional RBX1 binder. This negative endpoint is important: it reinforces the gap PRIB-KI is designed to address rather than weakening the project rationale.

## Workflow relevance

The case demonstrated an end-to-end AI-to-lab path:

```text
RBX1 target preparation
    -> de novo binder generation
    -> computational filtering and candidate comparison
    -> PRIB-style reliability and developability triage
    -> competition submission
    -> external selection
    -> wet-lab expression and binding evaluation
```

The design workflow used modern protein-design tooling, including **NVIDIA Proteina-Complexa** and **BindCraft**, followed by candidate filtering and reliability-oriented review. PRIB's role was not to be another generative model. Its contribution was the downstream question that later became central to PRIB-KI:

```text
Which computationally plausible protein designs deserve experimental resources first?
```

## Project relevance

This case supports the project in several specific ways:

- the team has operated a modern AI protein-design workflow under external competition constraints;
- the team has used state-of-the-art generative design tools and candidate-screening logic in practice;
- the submitted design set achieved strong external selection performance;
- selected candidates progressed into wet-lab testing rather than stopping at in-silico ranking;
- the team has direct experience with the gap between computationally attractive designs and experimental success;
- the negative binder outcome motivates a stronger reliability and developability layer before experimental allocation.

## Claim boundaries

This case must not be used to overclaim.

It does not prove:

- PRIB-KI has validated predictive accuracy;
- PRIB-KI has already enriched physical failures in a prospective reliability study;
- the selected designs were experimentally confirmed binders;
- the current public PRIB-KI demo is calibrated on wet-lab results;
- an external organization endorsed or validated PRIB-KI as a product.

Preferred wording:

> In the GEM x Adaptyv RBX1 Binder Design Competition linked to the ICLR 2026 GEM Workshop, the team used modern AI protein-design tools, including Proteina-Complexa and BindCraft, and applied reliability-oriented candidate triage. Seven of 21 submitted candidates advanced to wet-lab testing, corresponding to approximately the top 2.5% of more than 12,000 global submissions. No confirmed binder was ultimately obtained, which reinforced the need for stronger post-generation reliability screening.

## Ecosystem signal

The work also received positive attention and engagement from researchers in the NVIDIA generative-biology ecosystem, including public social-media engagement and follow-up discussion interest.

This should be described carefully. Unless a formal written basis is available, external materials should not state that NVIDIA endorsed, validated, sponsored, partnered with, or officially approved PRIB-KI.

## Link to PRIB-KI

The strategic lesson from this case is simple:

```text
generation can create plausible candidates;
selection determines which candidates deserve experiments;
wet-lab results expose the limits of computational plausibility.
```

PRIB-KI is built around that lesson. It is a post-generation physical reliability and developability layer for AI-designed proteins.
