from __future__ import annotations
import numpy as np
import pandas as pd
 



N_UNITS = 1000
T_STEPS = 20000
N_SEEDS = 4
M_GRID  = [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.30]
SUBS    = [1.0, 0.10]
 

BASE = dict(h=0.20, tau_rec=300.0, u=0.05)
OAT  = {
    "h":       [0.10, 0.20, 0.40],
    "tau_rec": [150.0, 300.0, 600.0],
    "u":       [0.03, 0.05, 0.08],
}
 

ONSET_THR, ONSET_QUIET, GROW_S, CEIL_FRAC = 3, 6, 4, 0.10
 
 
def simulate_vec(m: float, h: float, tau_rec: float, u: float,
                 n_seeds: int = N_SEEDS, t_steps: int = T_STEPS,
                 n_units: int = N_UNITS, seed0: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed0)
    A = np.zeros((n_seeds, t_steps), dtype=np.int64)
    r = np.ones(n_seeds)
    for t in range(t_steps - 1):
        rate = m * r * A[:, t] * (1.0 - A[:, t] / n_units) + h
        A[:, t + 1] = np.minimum(n_units, rng.poisson(np.clip(rate, 0, None)))
        r = r + (1.0 - r) / tau_rec - u * A[:, t + 1] / n_units
        r = np.clip(r, 0.0, 1.0)
    return A
 
 
def onset_ensemble_m(A: np.ndarray, f: float, rng: np.random.Generator,
                     n_units: int = N_UNITS) -> float:
    obs = A if f >= 1.0 else rng.binomial(A, f)
    onsets = [t for t in range(ONSET_QUIET, len(A) - GROW_S)
              if A[t] >= ONSET_THR and np.all(A[t - ONSET_QUIET:t] < ONSET_THR)]
    if len(onsets) < 30:
        return np.nan
    ens = np.array([obs[o:o + GROW_S + 1] for o in onsets], dtype=float).mean(0)
    ceil = CEIL_FRAC * n_units * f
    idx = np.where((ens > 0) & (ens < ceil))[0]
    idx = idx[idx >= 1]
    if len(idx) < 3:
        return np.nan
    return float(np.exp(np.polyfit(idx, np.log(ens[idx]), 1)[0]))
 
 
def crossing_point(m_true: np.ndarray, m_hat: np.ndarray) -> float:
    ok = np.isfinite(m_hat)
    mt, mh = m_true[ok], m_hat[ok]
    for i in range(len(mt) - 1):
        if mh[i] <= 1.0 < mh[i + 1]:
            frac = (1.0 - mh[i]) / (mh[i + 1] - mh[i])
            return float(mt[i] + frac * (mt[i + 1] - mt[i]))
    return np.nan
 
 
def m_onset_curve(world: dict, f: float) -> np.ndarray:
    curve = []
    for m in M_GRID:
        A = simulate_vec(m, **world, seed0=hash((m, f)) % 2**31)
        rng = np.random.default_rng(777)
        vals = [onset_ensemble_m(A[s], f, rng) for s in range(A.shape[0])]
        curve.append(np.nanmean(vals))
    return np.array(curve)
 
 
def main() -> None:
    mt = np.array(M_GRID)
    rows = []
    print(f"{'param':>9} {'value':>7} | {'m*(f=1)':>8} {'m*(f=.1)':>9} "
          f"{'monotone':>9} {'allsuper>1':>11}")
    print("-" * 60)
    for param, values in OAT.items():
        for v in values:
            world = dict(BASE); world[param] = v
            c1 = m_onset_curve(world, 1.0)
            c2 = m_onset_curve(world, 0.10)
            star1, star2 = crossing_point(mt, c1), crossing_point(mt, c2)
            mono = bool(np.all(np.diff(c1[np.isfinite(c1)]) > -0.02))
            allsuper = bool(np.all(c1[mt > 1.0] > 1.0))
            rows.append(dict(param=param, value=v, crossing_f1=star1,
                             crossing_f01=star2, monotone=mono, allsuper_gt1=allsuper,
                             **{f"m_{m:.2f}": c1[i] for i, m in enumerate(M_GRID)}))
            base_tag = "  <-base" if v == BASE[param] else ""
            print(f"{param:>9} {v:>7.2f} | {star1:>8.3f} {star2:>9.3f} "
                  f"{str(mono):>9} {str(allsuper):>11}{base_tag}")
 
    df = pd.DataFrame(rows)
    out = "transient_onset_paramsweep.csv"
    df.to_csv(out, index=False)
    print(f"\nSaved -> {out}")
 
    stars = df.crossing_f1.dropna().values
    print("\n--- READOUT ---")
    print(f"Crossing m* across all {len(stars)} worlds (full coverage): "
          f"{stars.min():.3f} to {stars.max():.3f}  (ideal = 1.000)")
    print(f"Max deviation of the crossing from 1.0: {np.max(np.abs(stars - 1)):.3f}")
    print(f"Monotone in every world: {bool(df.monotone.all())};  "
          f"all supercritical read >1 in every world: {bool(df.allsuper_gt1.all())}")
    drift = np.max(np.abs(stars - 1))
    if drift <= 0.05:
        print("=> Crossing is STABLE (<=0.05 drift). The sub/supercritical classifier "
              "is robust to the saturation/depression model.\n   The estimator is "
              "safe to take to real culture data (with the magnitude-compression "
              "caveat from step10 still standing).")
    else:
        print(f"=> Crossing DRIFTS by up to {drift:.3f}. Identify which parameter "
              "moves it (see table) -- that is the estimator's limit and must be\n"
              "   reported as a condition on the datasets where it can be trusted.")
 
 
if __name__ == "__main__":
    main()