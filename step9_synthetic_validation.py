"""
step9_synthetic_mr_validation.py
============================================================================
Synthetic supercritical validation harness for the criticality estimators,
BEFORE any in-vitro MEA data is touched.
 
WHY THIS EXISTS
---------------
The whole point of the slice/culture arm is to push the branching ratio m
above 1 so the MPR_C complexity peak finally gets bracketed. But the estimator
that would measure that -- Wilting-Priesemann MR -- has only ever been trusted
here in the m <= 1 regime, and per CLAUDE.md Section 2 every criticality
function needs a synthetic ground-truth check before it's used on real data.
There is currently NO such check for MR against known m > 1. This script is it.
 
It also asks the second question the EEG work raised: naive avalanche sigma was
garbage on EEG *because of severe subsampling*. Dense MEA is the least-subsampled
regime available -- so avalanche-based measures (sigma, and the avalanche-size
exponent tau) may be rehabilitated there and may see supercriticality that a
stationary autocorrelation fit cannot. This harness runs all three head to head.
 
MODEL (culture-like, saturating branching process)
--------------------------------------------------
Discrete-time population activity A(t) on N units with logistic saturation:
 
    A(t+1) ~ Poisson( m * A(t) * (1 - A(t)/N) + h )
 
- m       generative branching ratio (the ground truth we try to recover)
- N       population size (saturation ceiling)
- h       small external drive (keeps subcritical runs from going extinct)
- logistic (1 - A/N) damping = finite neurons + refractoriness; for m>1 the
  process rides at a STABLE fixed point A* = N(1 - 1/m) instead of exploding.
 
IMPORTANT, STATED UP FRONT (this is a prediction, verified below, not a finding
smuggled in): for a *stationary* saturated system the locally-linearized
branching ratio at the operating point is <= 1 even when generative m > 1.
For this logistic model the local slope is exactly (2 - m) on the supercritical
branch. So a stationary autocorrelation estimator like MR is expected to FOLD
the supercritical branch back below 1. If that happens, MR run on steady-state
culture activity cannot report m > 1, and the slice arm needs avalanche-shape
analysis of transients (or the unsaturated onset) instead. The harness's job is
to show cleanly whether that's true.
 
Subsampling model: observe each active unit independently with probability f
(binomial subsampling) -- the exact model MR was designed to be robust to.
 
Run:
    pip install --break-system-packages numpy scipy pandas mrestimator
    python3 step9_synthetic_mr_validation.py
Writes: synthetic_mr_validation.csv
============================================================================
"""
from __future__ import annotations
import logging
import numpy as np
import pandas as pd
 
logging.getLogger("mrestimator").setLevel(logging.ERROR)
 
# ---------------------------------------------------------------------------
# CONFIG (fixed before looking at results)
# ---------------------------------------------------------------------------
N_UNITS   = 1000
DRIVE_H   = 2.0
T_STEPS   = 8000
BURN_IN   = 2000
N_TRIALS  = 8
KMAX      = 150
M_GRID    = [0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.30]
SUBSAMPLE = [1.0, 0.10]          # full, and 10% MEA-like coverage
SEED      = 20260706
 
 
def simulate_branching(m: float, n_units: int = N_UNITS, t_steps: int = T_STEPS,
                       n_trials: int = N_TRIALS, h: float = DRIVE_H,
                       burn_in: int = BURN_IN, rng: np.random.Generator | None = None
                       ) -> np.ndarray:
    """Saturating (logistic) branching process. Returns (n_trials, t_steps-burn_in)
    integer activity, vectorized across trials (only time is sequential)."""
    rng = rng or np.random.default_rng()
    A = np.full(n_trials, max(1.0, h / max(1e-3, 1 - m)) if m < 1 else n_units * (1 - 1 / m))
    A = np.clip(A, 1, n_units).astype(float)
    out = np.empty((n_trials, t_steps), dtype=np.int64)
    for t in range(t_steps):
        rate = m * A * (1.0 - A / n_units) + h
        rate = np.clip(rate, 0, None)
        A = rng.poisson(rate).astype(float)
        A = np.clip(A, 0, n_units)
        out[:, t] = A
    return out[:, burn_in:]
 
 
def subsample(activity: np.ndarray, f: float, rng: np.random.Generator) -> np.ndarray:
    """Binomial subsampling: observe each active unit with probability f."""
    if f >= 1.0:
        return activity
    return rng.binomial(activity.astype(np.int64), f)
 
 
def mr_estimate(trials: np.ndarray, kmax: int = KMAX) -> tuple[float, float]:
    """Wilting-Priesemann MR (trialseparated), same call as step3b. -> (m_hat, tau)."""
    import mrestimator as mre
    if trials.sum() < 50:
        return np.nan, np.nan
    try:
        rk = mre.coefficients(trials, steps=(1, kmax), dt=1, method="trialseparated")
        ft = mre.fit(rk, fitfunc="exponential")
        return float(ft.mre), float(ft.tau)
    except Exception:
        return np.nan, np.nan
 
 
def avalanche_sigma_naive(trials: np.ndarray) -> float:
    """BIASED REFERENCE. Mean bin-to-bin ratio within avalanches (excursions above
    the per-trial median). Same spirit as step3's naive sigma."""
    sig = []
    for a in trials:
        thr = np.median(a)
        A = np.where(a > thr, a - thr, 0.0)
        b = 0
        while b < len(A):
            if A[b] > 0:
                st = b
                while b < len(A) and A[b] > 0:
                    b += 1
                av = A[st:b]
                if len(av) >= 2 and av[0] > 0:
                    sig.append(av[1] / av[0])
            else:
                b += 1
    return float(np.mean(sig)) if sig else np.nan
 
 
def avalanche_tau(trials: np.ndarray) -> tuple[float, int]:
    """Avalanche-size power-law exponent via continuous MLE (Clauset-Shalizi-Newman).
    Avalanche = contiguous run above per-trial median; size = summed excess activity.
    -> (tau, n_avalanches). Mean-field critical reference is tau ~ 1.5."""
    sizes = []
    for a in trials:
        thr = np.median(a)
        A = np.where(a > thr, a - thr, 0.0)
        b = 0
        while b < len(A):
            if A[b] > 0:
                st = b
                while b < len(A) and A[b] > 0:
                    b += 1
                sizes.append(A[st:b].sum())
            else:
                b += 1
    sizes = np.array([s for s in sizes if s >= 1.0])
    if len(sizes) < 30:
        return np.nan, len(sizes)
    xmin = 1.0
    tau = 1.0 + len(sizes) / np.sum(np.log(sizes / (xmin * 0.5)))
    return float(tau), int(len(sizes))
 
 
def local_slope_prediction(m: float) -> float:
    """Analytic local branching ratio at the model's operating point:
    subcritical branch (m<=1): ~ m ; supercritical branch (m>1): ~ 2 - m."""
    return m if m <= 1.0 else 2.0 - m
 
 
def main() -> None:
    rng = np.random.default_rng(SEED)
    rows = []
    print(f"{'m_true':>7} {'f':>5} {'m_MR':>7} {'pred':>6} {'tau_bins':>8} "
          f"{'sigma':>7} {'aval_tau':>8} {'n_aval':>7}")
    print("-" * 62)
    for m in M_GRID:
        full = simulate_branching(m, rng=rng)
        for f in SUBSAMPLE:
            obs = subsample(full, f, rng)
            m_hat, tau = mr_estimate(obs)
            sigma = avalanche_sigma_naive(obs)
            atau, n_av = avalanche_tau(obs)
            pred = local_slope_prediction(m)
            rows.append(dict(m_true=m, subsample_f=f, m_MR=m_hat, tau_bins=tau,
                             local_slope_pred=pred, sigma_naive=sigma,
                             aval_tau=atau, n_avalanches=n_av,
                             mean_activity=float(obs.mean())))
            print(f"{m:>7.2f} {f:>5.2f} {m_hat:>7.3f} {pred:>6.2f} "
                  f"{tau:>8.1f} {sigma:>7.3f} {atau:>8.3f} {n_av:>7d}")
 
    df = pd.DataFrame(rows)
    out = "synthetic_mr_validation.csv"
    df.to_csv(out, index=False)
    print(f"\nSaved -> {out}")
 
    # verdict
    sudf = df[(df.m_true > 1.0) & (df.subsample_f == 1.0)]
    max_super_mr = sudf.m_MR.max() if len(sudf) else np.nan
    print("\n--- READOUT ---")
    print(f"Max m_MR reported on genuinely supercritical (m_true>1) full-coverage "
          f"runs: {max_super_mr:.3f}")
    if np.isfinite(max_super_mr) and max_super_mr < 1.0:
        print("=> MR does NOT report m>1 on stationary saturated data. It folds the "
              "supercritical branch below 1, as predicted.\n   Implication: MR on "
              "steady-state culture activity cannot bracket the peak. The slice arm "
              "needs avalanche-shape / transient analysis, not stationary MR.")
    print("Check whether aval_tau moves monotonically across the m grid -- if it "
          "does, avalanche size distribution is the estimator that carries\n"
          "supercritical information on dense data.")
 
 
if __name__ == "__main__":
    main()