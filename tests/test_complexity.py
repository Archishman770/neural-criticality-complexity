"""Synthetic-signal regression tests for the complexity measures (CLAUDE.md Section 2).

Pins `lz_complexity` and `mpr_complexity` from step2_batch_extract.py to the reference
values recorded in CLAUDE.md: LZ ~ 0.012 / 1.037 / 1.025 and MPR_C ~ 0.211 / 0.028 /
0.484 for sine / noise / chaos. The generating signals were not recorded with those
numbers; the ones below were reconstructed on 2026-09-27 and reproduce all six values
to three decimals (n = 5000; sine = 100 cycles over linspace; noise = N(0,1) from
default_rng(0); chaos = logistic map r = 4, x0 = 0.4).

Run:  python3 -m pytest tests/
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import step2_batch_extract as s2  # noqa: E402

N = 5000
TOL = 0.002


def sine() -> np.ndarray:
    return np.sin(np.linspace(0, 2 * np.pi * 100, N))


def noise() -> np.ndarray:
    return np.random.default_rng(0).standard_normal(N)


def chaos(r: float = 4.0, x0: float = 0.4) -> np.ndarray:
    x = np.empty(N)
    x[0] = x0
    for i in range(1, N):
        x[i] = r * x[i - 1] * (1 - x[i - 1])
    return x


SIGNALS = {"sine": sine, "noise": noise, "chaos": chaos}
REF_LZ = {"sine": 0.012, "noise": 1.037, "chaos": 1.025}
REF_MPR = {"sine": 0.211, "noise": 0.028, "chaos": 0.484}


@pytest.mark.parametrize("name", SIGNALS)
def test_lz_reference(name: str) -> None:
    assert s2.lz_complexity(SIGNALS[name]()) == pytest.approx(REF_LZ[name], abs=TOL)


@pytest.mark.parametrize("name", SIGNALS)
def test_mpr_c_reference(name: str) -> None:
    _, c = s2.mpr_complexity(SIGNALS[name]())
    assert c == pytest.approx(REF_MPR[name], abs=TOL)


def test_lz_cannot_separate_noise_from_chaos() -> None:
    """LZ is a randomness measure: noise and chaos look alike (Section 2)."""
    assert abs(s2.lz_complexity(noise()) - s2.lz_complexity(chaos())) < 0.05


def test_mpr_c_separates_all_three() -> None:
    """MPR_C is non-monotonic in disorder: chaos > sine > noise (Section 2)."""
    c = {k: s2.mpr_complexity(f())[1] for k, f in SIGNALS.items()}
    assert c["chaos"] > c["sine"] > c["noise"]


def test_mpr_entropy_bounds() -> None:
    """Normalized permutation entropy H lies in [0, 1] and is near 1 for noise.

    Not exactly 1: d = 6 gives 720 patterns from ~5000 windows, so finite-sample
    bias pulls H down (0.989 here). MPR_C for noise (0.028) is that same bias.
    """
    h, _ = s2.mpr_complexity(noise())
    assert 0.0 <= h <= 1.0
    assert h > 0.98


def test_inline_lz76_matches_antropy() -> None:
    """The inline LZ76 fallback must equal antropy, so the result does not depend on
    whether antropy happens to be installed."""
    ant = pytest.importorskip("antropy")
    for f in SIGNALS.values():
        x = f()
        b = (x > np.median(x)).astype(int)
        inline = s2.lz76(b) * np.log2(len(b)) / len(b)
        assert inline == pytest.approx(ant.lziv_complexity(b, normalize=True), rel=1e-12)
