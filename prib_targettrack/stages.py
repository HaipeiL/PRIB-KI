"""Canonical, shared stages for the historical TargetTrack funnel."""

CANONICAL_STAGES = (
    "selected",
    "cloned",
    "expressed",
    "soluble",
    "purified",
    "crystallized",
    "diffraction_quality",
    "structure_determined",
    "deposited",
)

STAGE_RANK = {stage: rank for rank, stage in enumerate(CANONICAL_STAGES)}


def stage_at(rank: int) -> str:
    return CANONICAL_STAGES[rank]


def transition_name(rank: int) -> str:
    return f"{stage_at(rank)}_to_{stage_at(rank + 1)}"
