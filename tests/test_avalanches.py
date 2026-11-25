"""Ground-truth tests for src/criticality/avalanches.py.

A Galton-Watson branching process with Poisson(m) offspring has known mean-field
avalanche exponents at m = 1: sizes tau = 3/2, durations alpha = 2, and the
scaling relation beta = (alpha - 1)/(tau - 1) = 2 holds, so DCC ~ 0 and kappa ~ 1.
Subcritical m gives kappa < 1.

Run:  python3 -m pytest tests/
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.criticality.avalanches import dcc, extract_avalanches, fit_powerlaw, kappa  # noqa: E402

MAX_GEN = 2000
MAX_SIZE = 200_000


def galton_watson(m: float, n_aval: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Sizes and durations of n_aval single-seed Poisson(m) branching avalanches.

    Avalanches hitting MAX_GEN or MAX_SIZE are dropped (censored), as their true
    extent is unknown.
    """
    rng = np.random.default_rng(seed)
    sizes, durs = [], []
    for _ in range(n_aval):
        n, s, g = 1, 1, 1
        while n > 0 and g < MAX_GEN and s < MAX_SIZE:
            n = rng.poisson(m * n)
            s += n
            g += int(n > 0)
        if n == 0:
            sizes.append(s)
            durs.append(g)
    return np.array(sizes), np.array(durs)


def test_extract_avalanches() -> None:
    counts = np.array([3, 0, 1, 2, 0, 0, 5, 0, 4, 4, 4, 0, 7])
    s, d = extract_avalanches(counts)
    # leading [3] and trailing [7] touch the edges and are dropped
    assert s.tolist() == [3, 5, 12]
    assert d.tolist() == [2, 1, 3]


@pytest.mark.parametrize("a", [1.5, 2.5])
def test_fit_powerlaw_recovers_exponent(a: float) -> None:
    x = np.random.default_rng(1).zipf(a, 20_000)
    x = x[x <= 10_000]
    alpha, xmin, _ = fit_powerlaw(x, xmax=10_000)
    assert alpha == pytest.approx(a, abs=0.06)


def test_critical_galton_watson_exponents() -> None:
    s, d = galton_watson(1.0, 20_000, seed=2)
    r = dcc(s, d)
    assert r["tau"] == pytest.approx(1.5, abs=0.1)
    assert r["alpha"] == pytest.approx(2.0, abs=0.2)
    assert r["dcc"] < 0.25


def test_kappa_near_one_at_criticality_and_below_when_subcritical() -> None:
    s1, _ = galton_watson(1.0, 20_000, seed=3)
    s0, _ = galton_watson(0.8, 20_000, seed=4)
    k1, k0 = kappa(s1), kappa(s0)
    assert k1 == pytest.approx(1.0, abs=0.1)
    assert k0 < k1
    assert k0 < 1.0
