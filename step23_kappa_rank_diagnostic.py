"""Rank diagnostic for kappa across heterogeneous synthetic recordings (A9 §9.1).

Implements the specification fixed in docs/A9_culture_arm_preregistration.md §9.1
on 2026-10-03, before this script was first run. Do not change constants here
without logging an amendment.

Each synthetic recording draws its own branching ratio m and network parameters, and
observes a fixed random subset of units with lognormal propensities, so the effective
observation fraction p_rec differs across recordings (the main source of the step22
rank reversals). Reports how well kappa ranks m overall, within p_rec tertiles, near
the transition (|m - 1| <= 0.10), as a reversal rate, and as a partial correlation
controlling for p_rec and mean observed rate.

This is a diagnostic under the stated model assumptions, not a biological calibration.

Usage
-----
    python3 step23_kappa_rank_diagnostic.py
Writes kappa_rank_diagnostic_recordings.csv and kappa_rank_diagnostic_summary.csv.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from src.criticality.avalanches import extract_avalanches, kappa

SEED = 20261003
N_REC = 400
N_UNITS = 2000
T_STEPS = 300_000            # 1 step = 1 ms -> 5 min
BURN_IN = 5_000
M_RANGE = (0.80, 1.30)
H_RANGE = (0.10, 0.40)
TAU_REC_RANGE = (150.0, 600.0)
U_RANGE = (0.03, 0.08)
N_OBS_RANGE = (20, 142)      # hc-8 sorted-neuron range
PROPENSITY_SIGMA = 1.0       # w_i ~ LogNormal(0, 1)
BINS_MS = [2, 4, 8]
PRIMARY_BIN = 4
MIN_AVAL = 500
NEAR = 0.10                  # |m - 1| <= NEAR
DM_PAIR = 0.10               # pairs with |dm| >= DM_PAIR for the reversal rate
CRITERIA = dict(rho_overall=0.80, rho_tertile=0.70, rho_near=0.50, reversal_max=0.10)


def draw_population(rng: np.random.Generator) -> pd.DataFrame:
    """Per-recording parameters and effective observation fraction p_rec.

    Args:
        rng: random generator.

    Returns:
        DataFrame with m, h, tau_rec, u, n_obs, p_rec per recording.
    """
    w = rng.lognormal(0.0, PROPENSITY_SIGMA, N_UNITS)
    n_obs = rng.integers(N_OBS_RANGE[0], N_OBS_RANGE[1] + 1, N_REC)
    p_rec = np.array([w[rng.choice(N_UNITS, k, replace=False)].sum() / w.sum() for k in n_obs])
    return pd.DataFrame(dict(
        rec=np.arange(N_REC),
        m=rng.uniform(*M_RANGE, N_REC),
        h=rng.uniform(*H_RANGE, N_REC),
        tau_rec=rng.uniform(*TAU_REC_RANGE, N_REC),
        u=rng.uniform(*U_RANGE, N_REC),
        n_obs=n_obs,
        p_rec=p_rec,
    ))


def simulate_observed(pop: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    """Run all recordings in parallel (vectorized over recordings; time is sequential).

    Same update as step11.simulate_vec, with per-recording parameter arrays. The
    observed count is Binomial(A(t), p_rec), the exact marginal of spikes assigned to
    units in proportion to fixed propensities.

    Returns:
        (N_REC, T_STEPS) uint16 observed counts after burn-in.
    """
    m, h = pop.m.values, pop.h.values
    tau, u, p = pop.tau_rec.values, pop.u.values, pop.p_rec.values
    A = np.zeros(N_REC)
    r = np.ones(N_REC)
    obs = np.empty((N_REC, T_STEPS), dtype=np.uint16)
    for t in range(BURN_IN + T_STEPS):
        rate = m * r * A * (1.0 - A / N_UNITS) + h
        A = np.minimum(N_UNITS, rng.poisson(np.clip(rate, 0.0, None))).astype(float)
        r = np.clip(r + (1.0 - r) / tau - u * A / N_UNITS, 0.0, 1.0)
        if t >= BURN_IN:
            obs[:, t - BURN_IN] = rng.binomial(A.astype(np.int64), p)
    return obs


def rebin(x: np.ndarray, k: int) -> np.ndarray:
    """Sum consecutive groups of k bins (drops the remainder)."""
    n = (len(x) // k) * k
    return x[:n].reshape(-1, k).sum(1) if k > 1 else x


def reversal_rate(m: np.ndarray, k: np.ndarray, dm: float) -> float:
    """Fraction of pairs with m_i - m_j >= dm where the higher-m recording has lower kappa."""
    dmat = m[:, None] - m[None, :]
    dk = k[:, None] - k[None, :]
    sel = dmat >= dm                       # i has the higher m
    return float(np.mean(dk[sel] < 0)) if sel.any() else np.nan


def partial_spearman(x: np.ndarray, y: np.ndarray, covs: np.ndarray) -> float:
    """Spearman correlation of x and y after regressing both ranks on covariate ranks."""
    rx, ry = rankdata(x), rankdata(y)
    C = np.column_stack([np.ones(len(x))] + [rankdata(c) for c in covs.T])
    ex = rx - C @ np.linalg.lstsq(C, rx, rcond=None)[0]
    ey = ry - C @ np.linalg.lstsq(C, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def main() -> None:
    """Simulate, compute kappa, save per-recording table, then summarize.

    With `--summary-only`, skip simulation and summarize the saved recordings CSV.
    """
    import sys
    if "--summary-only" in sys.argv:
        summarize(pd.read_csv("kappa_rank_diagnostic_recordings.csv"))
        return
    rng = np.random.default_rng(SEED)
    pop = draw_population(rng)
    print(f"simulating {N_REC} recordings x {T_STEPS} steps ...", flush=True)
    obs = simulate_observed(pop, rng)
    pop["mean_obs_rate"] = obs.mean(1)

    for b in BINS_MS:
        ks, na = [], []
        for i in range(N_REC):
            sz, _ = extract_avalanches(rebin(obs[i].astype(np.int64), b))
            na.append(len(sz))
            ks.append(kappa(sz) if len(sz) >= MIN_AVAL else np.nan)
        pop[f"kappa_{b}"] = ks
        pop[f"n_aval_{b}"] = na
    pop.to_csv("kappa_rank_diagnostic_recordings.csv", index=False)
    summarize(pop)


def summarize(pop: pd.DataFrame) -> None:
    """Rank metrics per bin, the §9.1 criteria at the primary bin, and the summary CSV.

    Args:
        pop: per-recording table with m, p_rec, mean_obs_rate and kappa_{bin} columns.
    """
    tert = pd.qcut(pop.p_rec, 3, labels=["low", "mid", "high"])
    rows = []
    for b in BINS_MS:
        d = pop.assign(tert=tert).dropna(subset=[f"kappa_{b}"])
        k, m = d[f"kappa_{b}"].values, d.m.values
        near = np.abs(m - 1.0) <= NEAR
        row = dict(bin_ms=b, n_used=len(d), n_excluded=N_REC - len(d),
                   rho_overall=spearmanr(k, m)[0],
                   rho_near=spearmanr(k[near], m[near])[0] if near.sum() > 5 else np.nan,
                   n_near=int(near.sum()),
                   reversal_overall=reversal_rate(m, k, DM_PAIR),
                   reversal_near=reversal_rate(m[near], k[near], DM_PAIR),
                   partial_rho=partial_spearman(k, m, d[["p_rec", "mean_obs_rate"]].values))
        for t in ["low", "mid", "high"]:
            s = d[d.tert == t]
            row[f"rho_tertile_{t}"] = spearmanr(s[f"kappa_{b}"], s.m)[0] if len(s) > 5 else np.nan
        rows.append(row)
    summ = pd.DataFrame(rows)
    summ.to_csv("kappa_rank_diagnostic_summary.csv", index=False)

    print("\n=== kappa rank diagnostic (A9 §9.1) ===")
    print(summ.round(3).T.to_string(header=False))
    p = summ[summ.bin_ms == PRIMARY_BIN].iloc[0]
    checks = {
        f"overall rho >= {CRITERIA['rho_overall']}": p.rho_overall >= CRITERIA["rho_overall"],
        f"every p_rec tertile rho >= {CRITERIA['rho_tertile']}":
            min(p.rho_tertile_low, p.rho_tertile_mid, p.rho_tertile_high) >= CRITERIA["rho_tertile"],
        f"near-transition rho >= {CRITERIA['rho_near']}": p.rho_near >= CRITERIA["rho_near"],
        f"overall reversal <= {CRITERIA['reversal_max']}": p.reversal_overall <= CRITERIA["reversal_max"],
    }
    print(f"\n--- criteria at {PRIMARY_BIN} ms ---")
    for name, ok in checks.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    verdict = all(checks.values())
    print(f"\nDIAGNOSTIC {'PASSED' if verdict else 'FAILED'} (under the stated model assumptions; "
          "kappa remains an avalanche-distribution axis either way)")
    print(f"p_rec range {pop.p_rec.min():.4f}-{pop.p_rec.max():.4f}; "
          f"tertile cut points {np.quantile(pop.p_rec, [1/3, 2/3]).round(4)}")


if __name__ == "__main__":
    main()
