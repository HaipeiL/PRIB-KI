"""Small dependency-free intervals for reproducible progression summaries."""

from __future__ import annotations

import math


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return (math.nan, math.nan)
    p = successes / total
    denominator = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    radius = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total) / denominator
    return (max(0.0, centre - radius), min(1.0, centre + radius))


def _beta_fraction(a: float, b: float, x: float) -> float:
    tiny = 3.0e-30
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = 1.0 / max(d, tiny)
    value = d
    for step in range(1, 201):
        twice = 2 * step
        aa = step * (b - step) * x / ((qam + twice) * (a + twice))
        d = 1.0 + aa * d
        d = 1.0 / max(d, tiny)
        c = max(1.0 + aa / c, tiny)
        value *= d * c
        aa = -(a + step) * (qab + step) * x / ((a + twice) * (qap + twice))
        d = 1.0 + aa * d
        d = 1.0 / max(d, tiny)
        c = max(1.0 + aa / c, tiny)
        delta = d * c
        value *= delta
        if abs(delta - 1.0) < 3.0e-12:
            break
    return value


def regularized_beta(x: float, a: float, b: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    prefix = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return prefix * _beta_fraction(a, b, x) / a
    return 1.0 - prefix * _beta_fraction(b, a, 1.0 - x) / b


def beta_quantile(probability: float, a: float, b: float) -> float:
    lower, upper = 0.0, 1.0
    for _ in range(70):
        midpoint = (lower + upper) / 2.0
        if regularized_beta(midpoint, a, b) < probability:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def jeffreys_interval(successes: int, failures: int) -> tuple[float, float, float]:
    a, b = successes + 0.5, failures + 0.5
    return (a / (a + b), beta_quantile(0.025, a, b), beta_quantile(0.975, a, b))
