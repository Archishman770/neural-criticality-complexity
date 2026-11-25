"""Validate kappa and DCC on the bursting model under observed-only conditions (pre-A9).

The culture-arm criticality axis (user decision 2026-10-02) is Shew's kappa (signed:
<1 subcritical, >1 supercritical) plus DCC (unsigned distance from the crackling-noise
scaling relation). Both are known to shift under subsampling and with bin width.
Before any hc-8 / Charlesworth data is inspected, this script measures, on the same
saturating, depressing branching model used in step10/11/21:
  - where kappa crosses 1 as a function of true m (ideal: m = 1),
  - where DCC is minimal (ideal: m = 1),
  - how both move with observation fraction f and with bin width (in model steps).
Its output fixes what the culture arm may claim, and is recorded in the A9 plan.

Usage
-----
    python3 step22_avalanche_observed_validation.py
Writes avalanche_observed_validation.csv and avalanche_observed_summary.csv.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from src.criticality.avalanches import dcc, extract_avalanches, kappa
from step11_param_robustness_sweep import crossing_point, simulate_vec
from step21_transient_observed_validation import worlds

T_STEPS = 100_000
M_GRID = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.30]
F_OBS = [1.0, 0.10, 0.05, 0.02]
BIN_STEPS = [1, 2, 4]           # bin width in model steps (real bins != generation time)
MIN_AVAL = 200


def rebin(x: np.ndarray, k: int) -> np.ndarray:
    """Sum consecutive groups of k bins (drops the remainder)."""
    n = (len(x) // k) * k
    return x[:n].reshape(-1, k).sum(1) if k > 1 else x


def main() -> None:
    rows = []
    for wname, w in worlds():
        for m in M_GRID:
            A_all = simulate_vec(m, **w, t_steps=T_STEPS, seed0=int(1000 * m) + 7)
            for f in F_OBS:
                rng = np.random.default_rng(int(1e4 * f) + int(100 * m) + 3)
                obs_all = A_all if f >= 1.0 else rng.binomial(A_all, f)
                for k in BIN_STEPS:
                    ks, all_s, all_d = [], [], []
                    for s in range(obs_all.shape[0]):
                        sz, du = extract_avalanches(rebin(obs_all[s], k))
                        if len(sz) >= MIN_AVAL:
                            ks.append(kappa(sz))
                        all_s.append(sz)
                        all_d.append(du)
                    sz, du = np.concatenate(all_s), np.concatenate(all_d)
                    d = dcc(sz, du) if len(sz) >= MIN_AVAL else dict(
                        tau=np.nan, alpha=np.nan, beta_fit=np.nan, beta_pred=np.nan, dcc=np.nan)
                    rows.append(dict(world=wname, f=f, bin=k, m_true=m,
                                     kappa=np.nanmean(ks) if ks else np.nan,
                                     kappa_sd=np.nanstd(ks) if ks else np.nan,
                                     n_aval=len(sz), **d))
        print(f"done world {wname}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv("avalanche_observed_validation.csv", index=False)

    mt = np.array(M_GRID)
    summ = []
    for (wname, f, k), g in df.groupby(["world", "f", "bin"]):
        g = g.sort_values("m_true")
        kap = g.kappa.values
        dc = g.dcc.values
        ok = np.isfinite(kap) & (mt != 1.0)
        summ.append(dict(
            world=wname, f=f, bin=k,
            kappa_cross=crossing_point(mt, kap),
            kappa_sign_acc=(float(np.mean(np.sign(kap[ok] - 1) == np.sign(mt[ok] - 1)))
                            if ok.any() else np.nan),
            kappa_monotone=bool(np.all(np.diff(kap[np.isfinite(kap)]) > 0)),   # strict
            kappa_n_decreases=int(np.sum(np.diff(kap[np.isfinite(kap)]) <= 0)),
            dcc_argmin_m=float(mt[np.nanargmin(dc)]) if np.isfinite(dc).any() else np.nan,
            min_aval=int(g.n_aval.min())))
    sm = pd.DataFrame(summ)
    sm.to_csv("avalanche_observed_summary.csv", index=False)

    print("\n=== kappa crossing (ideal m = 1.00), sign accuracy, DCC argmin, by f and bin ===")
    agg = sm.groupby(["f", "bin"]).agg(
        cross_min=("kappa_cross", "min"), cross_max=("kappa_cross", "max"),
        cross_nan=("kappa_cross", lambda s: int(s.isna().sum())),
        sign_acc_min=("kappa_sign_acc", "min"), monotone_all=("kappa_monotone", "all"),
        dcc_argmin_lo=("dcc_argmin_m", "min"), dcc_argmin_hi=("dcc_argmin_m", "max"),
        min_aval=("min_aval", "min"))
    print(agg.round(3).to_string())
    print("\nbase world, kappa vs m_true (rows f, cols m):")
    b = df[(df.world == "base") & (df.bin == 1)].pivot_table(index="f", columns="m_true",
                                                              values="kappa")
    print(b.round(3).to_string())


if __name__ == "__main__":
    main()
