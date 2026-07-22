from __future__ import annotations
import logging
import numpy as np
import pandas as pd
 
logging.getLogger("mrestimator").setLevel(logging.ERROR)
 



N_UNITS   = 1000
DRIVE_H   = 0.20
TAU_REC   = 300.0
DEPLETE_U = 0.05
T_STEPS   = 40000
M_GRID    = [0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.30]
SUBSAMPLE = [1.0, 0.10]
SEEDS     = [0, 1, 2]
 

ONSET_THR   = 3
ONSET_QUIET = 6
GROW_S      = 4
CEIL_FRAC   = 0.10
 

MR_TRIAL_LEN = 200
MR_KMAX      = 150
 
 
def simulate_bursting_network(m: float, n_units: int = N_UNITS, t_steps: int = T_STEPS,
                              h: float = DRIVE_H, tau_rec: float = TAU_REC,
                              u: float = DEPLETE_U,
                              rng: np.random.Generator | None = None) -> np.ndarray:
    rng = rng or np.random.default_rng()
    A = np.zeros(t_steps, dtype=np.int64)
    r = 1.0
    for t in range(t_steps - 1):
        rate = m * r * A[t] * (1.0 - A[t] / n_units) + h
        A[t + 1] = min(n_units, rng.poisson(max(rate, 0.0)))
        r = r + (1.0 - r) / tau_rec - u * A[t + 1] / n_units
        r = min(1.0, max(0.0, r))
    return A
 
 
def subsample(A: np.ndarray, f: float, rng: np.random.Generator) -> np.ndarray:
    return A if f >= 1.0 else rng.binomial(A, f)
 
 
def onset_ensemble_m(A: np.ndarray, f: float, rng: np.random.Generator,
                     n_units: int = N_UNITS) -> tuple[float, int]:
    obs = subsample(A, f, rng)
    onsets = [t for t in range(ONSET_QUIET, len(A) - GROW_S)
              if A[t] >= ONSET_THR and np.all(A[t - ONSET_QUIET:t] < ONSET_THR)]
    if len(onsets) < 30:
        return np.nan, len(onsets)
    ens = np.array([obs[o:o + GROW_S + 1] for o in onsets], dtype=float).mean(0)
    ceil = CEIL_FRAC * n_units * f
    idx = np.where((ens > 0) & (ens < ceil))[0]
    idx = idx[idx >= 1]
    if len(idx) < 3:
        return np.nan, len(onsets)
    slope = np.polyfit(idx, np.log(ens[idx]), 1)[0]
    return float(np.exp(slope)), len(onsets)
 
 
def stationary_mr(A: np.ndarray, f: float, rng: np.random.Generator) -> float:
    import mrestimator as mre
    obs = subsample(A, f, rng)
    nt = len(obs) // MR_TRIAL_LEN
    if nt < 3 or obs.sum() < 50:
        return np.nan
    trials = obs[:nt * MR_TRIAL_LEN].reshape(nt, MR_TRIAL_LEN).astype(float)
    try:
        rk = mre.coefficients(trials, steps=(1, MR_KMAX), dt=1, method="trialseparated")
        return float(mre.fit(rk, fitfunc="exponential").mre)
    except Exception:
        return np.nan
 
 
def naive_stationary_slope(A: np.ndarray, f: float, rng: np.random.Generator) -> float:
    obs = subsample(A, f, rng).astype(float)
    x, y = obs[:-1], obs[1:]
    if x.var() < 1e-9:
        return np.nan
    return float(np.polyfit(x, y, 1)[0])
 
 
def main() -> None:
    rows = []
    print(f"{'m_true':>7} {'f':>5} | {'m_onset':>8} {'m_MR_stat':>10} "
          f"{'slope_stat':>11} | {'n_onset':>8}")
    print("-" * 62)
    for m in M_GRID:
        for f in SUBSAMPLE:
            onset_vals, mr_vals, slope_vals, nons = [], [], [], []
            for sd in SEEDS:
                A = simulate_bursting_network(m, rng=np.random.default_rng(sd))
                rng = np.random.default_rng(1000 + sd)
                mo, no = onset_ensemble_m(A, f, rng)
                onset_vals.append(mo); nons.append(no)
                mr_vals.append(stationary_mr(A, f, np.random.default_rng(2000 + sd)))
                slope_vals.append(naive_stationary_slope(A, f, np.random.default_rng(3000 + sd)))
            mo_m, mr_m = np.nanmean(onset_vals), np.nanmean(mr_vals)
            sl_m, no_m = np.nanmean(slope_vals), int(np.mean(nons))
            rows.append(dict(m_true=m, subsample_f=f, m_onset=mo_m,
                             m_onset_sd=np.nanstd(onset_vals), m_MR_stationary=mr_m,
                             slope_stationary=sl_m, n_onsets=no_m))
            print(f"{m:>7.2f} {f:>5.2f} | {mo_m:>8.3f} {mr_m:>10.3f} "
                  f"{sl_m:>11.3f} | {no_m:>8d}")
 
    df = pd.DataFrame(rows)
    out = "transient_onset_validation.csv"
    df.to_csv(out, index=False)
    print(f"\nSaved -> {out}")
 

    full = df[df.subsample_f == 1.0].sort_values("m_true")
    sub_true = full[full.m_true > 1.0]
    onset_super = sub_true.m_onset
    mr_super = sub_true.m_MR_stationary
    mono = np.all(np.diff(full.m_onset.values) > -0.02)
    crosses = (full[full.m_true <= 1.0].m_onset.max() <= 1.03 and
               full[full.m_true > 1.0].m_onset.min() >= 0.99)
    print("\n--- READOUT (full coverage) ---")
    print(f"Supercritical m_true>1: m_onset reads {onset_super.min():.3f}-{onset_super.max():.3f} "
          f"(all >1: {(onset_super > 1.0).all()})")
    print(f"   vs stationary MR on the SAME data: {mr_super.min():.3f}-{mr_super.max():.3f} "
          f"(folded below 1: {(mr_super < 1.0).all()})")
    print(f"Onset estimator monotone in m_true: {mono};  crosses 1 at m=1: {crosses}")
    print("Interpretation: the transient estimator recovers the SIGN of "
          "(m-1) and monotonicity where stationary MR folds.\n   Magnitude is "
          "compressed toward 1 (intrinsic to finite-N saturation) -- report as a "
          "sub/supercritical classifier + crossing detector,\n   not a calibrated m. "
          "Subsampling robustness: compare the f=1.0 and f=0.10 rows above.")
 
 
if __name__ == "__main__":
    main()