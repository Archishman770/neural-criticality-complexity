"""Avalanche-based criticality measures: extraction, power-law MLE, kappa, DCC.

Used for the culture arm (amendment A9), where stationary MR folds supercritical
values below 1 and the transient-onset estimator fails under realistic observation
(step21). All functions take a binned population count series; the bin width is a
caller decision and must be fixed in physical time across recordings (A8 lesson).

References
----------
Shew et al. 2009, J Neurosci 29:15595 (kappa).
Ma, Turrigiano, Wessel, Hengen 2019, Neuron 104:655; Fontenele et al. 2019 PRL
122:208101 (deviation from criticality coefficient, DCC; scaling relation
beta = (alpha - 1) / (tau - 1)).
Clauset, Shalizi, Newman 2009, SIAM Rev 51:661 (discrete power-law MLE, KS xmin).
"""
from __future__ import annotations

import numpy as np

ALPHA_GRID = np.linspace(1.01, 4.0, 600)


def extract_avalanches(counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split a binned count series into avalanches (runs of non-empty bins).

    An avalanche is a maximal run of consecutive bins with count > 0, bounded by
    empty bins. Runs touching either end of the series are dropped, since their
    true extent is unknown.

    Args:
        counts: 1-D non-negative integer counts per bin.

    Returns:
        (sizes, durations): total events and number of bins per avalanche.
    """
    counts = np.asarray(counts)
    active = counts > 0
    edges = np.diff(np.concatenate([[0], active.astype(np.int8), [0]]))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1)          # exclusive
    keep = (starts > 0) & (ends < len(active))
    starts, ends = starts[keep], ends[keep]
    csum = np.concatenate([[0], np.cumsum(counts)])
    return csum[ends] - csum[starts], ends - starts


def _loglik_grid(x: np.ndarray, xmin: int, xmax: int) -> np.ndarray:
    """Log-likelihood of a discrete power law truncated to [xmin, xmax], per alpha."""
    k = np.arange(xmin, xmax + 1, dtype=float)
    log_z = np.log(np.exp(-np.outer(ALPHA_GRID, np.log(k))).sum(1))
    return -ALPHA_GRID * np.log(x).sum() - len(x) * log_z


def fit_powerlaw(x: np.ndarray, xmax: int | None = None, xmin_candidates: int = 20,
                 min_tail: int = 50) -> tuple[float, int, float]:
    """Discrete truncated power-law MLE with xmin chosen by minimum KS distance.

    Args:
        x: positive integer samples (avalanche sizes or durations).
        xmax: upper truncation; defaults to max(x).
        xmin_candidates: number of distinct candidate xmin values tried (smallest first).
        min_tail: minimum samples in [xmin, xmax] for a candidate to be considered.

    Returns:
        (alpha, xmin, ks): exponent, chosen lower cutoff, KS distance of the fit.
        NaNs if no candidate has min_tail samples.
    """
    x = np.asarray(x, dtype=np.int64)
    x = x[x >= 1]
    xmax = int(xmax if xmax is not None else (x.max() if len(x) else 1))
    best = (np.nan, -1, np.inf)
    for xmin in np.unique(x)[:xmin_candidates]:
        tail = x[(x >= xmin) & (x <= xmax)]
        if len(tail) < min_tail or xmin >= xmax:
            continue
        a = float(ALPHA_GRID[np.argmax(_loglik_grid(tail, int(xmin), xmax))])
        k = np.arange(xmin, xmax + 1)
        pmf = k.astype(float) ** -a
        model_cdf = np.cumsum(pmf / pmf.sum())
        emp_cdf = np.searchsorted(np.sort(tail), k, side="right") / len(tail)
        ks = float(np.max(np.abs(emp_cdf - model_cdf)))
        if ks < best[2]:
            best = (a, int(xmin), ks)
    return best


def kappa(sizes: np.ndarray, tau_ref: float = 1.5, n_points: int = 10) -> float:
    """Shew et al. (2009) kappa: deviation of the size CDF from a reference power law.

    kappa = 1 + mean_k [F_ref(b_k) - F_data(b_k)] over n_points log-spaced sizes b_k
    between the smallest and largest avalanche. kappa > 1: excess large avalanches
    (supercritical side); kappa < 1: deficit (subcritical side).

    Args:
        sizes: avalanche sizes.
        tau_ref: exponent of the reference (critical) size distribution.
        n_points: number of log-spaced comparison points.

    Returns:
        kappa, or NaN if fewer than 2 distinct sizes.
    """
    s = np.asarray(sizes, dtype=float)
    s = s[s >= 1]
    if len(np.unique(s)) < 2:
        return np.nan
    lo, hi = s.min(), s.max()
    b = np.logspace(np.log10(lo), np.log10(hi), n_points)
    k = np.arange(int(lo), int(hi) + 1, dtype=float)
    pmf = k ** -tau_ref
    ref_cdf = np.cumsum(pmf / pmf.sum())
    f_ref = ref_cdf[np.clip(np.floor(b).astype(int) - int(lo), 0, len(k) - 1)]
    f_data = np.searchsorted(np.sort(s), b, side="right") / len(s)
    return float(1.0 + np.mean(f_ref - f_data))


def dcc(sizes: np.ndarray, durations: np.ndarray, min_tail: int = 50) -> dict:
    """Deviation from criticality coefficient, |beta_pred - beta_fit|.

    tau (sizes) and alpha (durations) by truncated discrete MLE; beta_fit is the
    slope of log mean size against log duration over the duration fit range;
    beta_pred = (alpha - 1) / (tau - 1) is the crackling-noise scaling relation.

    Args:
        sizes: avalanche sizes.
        durations: avalanche durations (bins), same length as sizes.
        min_tail: minimum samples for each power-law fit.

    Returns:
        dict with tau, alpha, beta_fit, beta_pred, dcc (NaNs if a fit fails).
    """
    s = np.asarray(sizes)
    d = np.asarray(durations)
    tau, _, _ = fit_powerlaw(s, min_tail=min_tail)
    alpha, dmin, _ = fit_powerlaw(d, min_tail=min_tail)
    out = dict(tau=tau, alpha=alpha, beta_fit=np.nan, beta_pred=np.nan, dcc=np.nan)
    if not (np.isfinite(tau) and np.isfinite(alpha)) or tau <= 1.0:
        return out
    sel = d >= dmin
    ud = np.unique(d[sel])
    if len(ud) < 3:
        return out
    mean_s = np.array([s[sel][d[sel] == u].mean() for u in ud])
    beta_fit = float(np.polyfit(np.log(ud), np.log(mean_s), 1)[0])
    beta_pred = (alpha - 1.0) / (tau - 1.0)
    out.update(beta_fit=beta_fit, beta_pred=beta_pred, dcc=abs(beta_pred - beta_fit))
    return out
