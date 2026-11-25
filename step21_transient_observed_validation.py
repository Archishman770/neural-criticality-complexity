"""Re-validate the transient-onset estimator under realistic observation (pre-A9).

step10/step11 detect burst onsets on the FULL population A (ground truth) and only
then average the subsampled trace. Real recordings never see A: in CRCNS hc-8 only
~91 sorted neurons out of thousands are observed. Detecting onsets on the observed
counts conditions on threshold crossings of a noisy trace, which can make growth
appear where there is none (selection / regression-to-the-mean). The fit ceiling in
step10/11 also uses N*f, and N is unknown in a culture.

This script runs both versions on identical simulations:
  - "oracle":   onsets from A, ceiling from N*f            (step10/11 behaviour)
  - "observed": onsets from obs, ceiling from observed burst peaks (deployable)
  - "split":    observed neurons split into two disjoint halves; onsets detected on
                one half, growth measured on the other, so selection on the
                detection half's noise cannot bias the measured growth (deployable)
across the step11 parameter worlds, three observation fractions, and two onset
thresholds. It also reports the growth rate lam = d ln<A>/dt per step, whose sign
(lam = 0 at criticality) does not depend on bin width, unlike m = exp(lam * dt).

Its output decides what the culture arm may claim (sign of m - 1 only, or offset
magnitude), and is to be fixed BEFORE any hc-8 / Charlesworth data is inspected.

Usage
-----
    python3 step21_transient_observed_validation.py
Writes transient_observed_validation.csv and transient_observed_crossings.csv.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from step11_param_robustness_sweep import BASE, OAT, crossing_point, simulate_vec

N_UNITS = 1000
T_STEPS = 100_000               # 5x step11: closer to the burst counts of a ~1 h recording
M_GRID = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.30]
F_OBS = [0.10, 0.05, 0.02]      # observed fraction; hc-8 is plausibly ~0.01-0.1
THR_OBS = [2, 3]                # onset threshold on the OBSERVED count
ONSET_QUIET = 6                 # bins below threshold required before an onset
GROW_S = 4                      # bins after onset used for the growth fit
CEIL_FRAC = 0.10                # fit only while ensemble mean < 10% of burst peak
PEAK_WIN = 50                   # bins after onset searched for the burst peak
MIN_ONSETS = 30
ORACLE_THR = 3                  # step10/11 onset threshold on A


def onsets_from(trace: np.ndarray, thr: int) -> np.ndarray:
    """Bins where trace >= thr after ONSET_QUIET bins below thr (vectorized).

    Args:
        trace: 1-D count series used for detection.
        thr: detection threshold in counts.

    Returns:
        Array of onset indices leaving room for the growth and peak windows.
    """
    above = trace >= thr
    below = ~above
    # quiet[t] is True if the ONSET_QUIET bins before t are all below threshold.
    c = np.concatenate([[0], np.cumsum(below)])
    t = np.arange(ONSET_QUIET, len(trace) - max(GROW_S, PEAK_WIN))
    quiet = (c[t] - c[t - ONSET_QUIET]) == ONSET_QUIET
    return t[above[t] & quiet]


def growth_fit(obs: np.ndarray, onsets: np.ndarray, ceil: float) -> tuple[float, float]:
    """Fit ln of the onset-aligned ensemble mean of obs over the unsaturated window.

    Args:
        obs: observed count series.
        onsets: onset indices.
        ceil: ensemble values at or above this are excluded (saturation).

    Returns:
        (lam per bin, m = exp(lam)); NaN if fewer than 3 usable points.
    """
    if len(onsets) < MIN_ONSETS:
        return np.nan, np.nan
    win = onsets[:, None] + np.arange(GROW_S + 1)[None, :]
    ens = obs[win].astype(float).mean(0)
    idx = np.where((ens > 0) & (ens < ceil))[0]
    idx = idx[idx >= 1]          # bin 0 is selected by the threshold itself
    if len(idx) < 3:
        return np.nan, np.nan
    lam = float(np.polyfit(idx, np.log(ens[idx]), 1)[0])
    return lam, float(np.exp(lam))


def estimate(A: np.ndarray, obs: np.ndarray, f: float, method: str, thr: int
             ) -> tuple[float, float, int]:
    """Transient-onset estimate under the oracle or the observed-only protocol.

    Args:
        A: full population activity (used only by the oracle protocol).
        obs: binomially subsampled activity.
        f: observation fraction (used only by the oracle ceiling).
        method: "oracle" or "observed".
        thr: onset threshold on obs (observed protocol only).

    Returns:
        (m_hat, lam, n_onsets)
    """
    if method == "oracle":
        on = onsets_from(A, ORACLE_THR)
        ceil = CEIL_FRAC * N_UNITS * f
    elif method == "split":
        # obs is a (2, T) pair of disjoint observed subsets: detect on [0], measure on [1].
        det, meas = obs[0], obs[1]
        on = onsets_from(det, thr)
        if len(on) < MIN_ONSETS:
            return np.nan, np.nan, len(on)
        peaks = meas[on[:, None] + np.arange(PEAK_WIN)[None, :]].max(1)
        ceil = max(CEIL_FRAC * float(np.median(peaks)), 1.0)
        lam, m = growth_fit(meas, on, ceil)
        return m, lam, len(on)
    else:
        on = onsets_from(obs, thr)
        if len(on) < MIN_ONSETS:
            return np.nan, np.nan, len(on)
        peaks = obs[on[:, None] + np.arange(PEAK_WIN)[None, :]].max(1)
        ceil = CEIL_FRAC * float(np.median(peaks))
        ceil = max(ceil, float(thr) + 1.0)   # never below the detection level
    lam, m = growth_fit(obs, on, ceil)
    return m, lam, len(on)


def worlds() -> list[tuple[str, dict]]:
    """Unique parameter worlds from step11's one-at-a-time sweep."""
    out = [("base", dict(BASE))]
    for p, vals in OAT.items():
        for v in vals:
            if v != BASE[p]:
                w = dict(BASE)
                w[p] = v
                out.append((f"{p}={v:g}", w))
    return out


def main() -> None:
    rows = []
    for wname, w in worlds():
        for m in M_GRID:
            A_all = simulate_vec(m, **w, t_steps=T_STEPS, seed0=int(1000 * m) + len(rows))
            for f in F_OBS:
                rng = np.random.default_rng(int(1e4 * f) + int(100 * m))
                obs_all = rng.binomial(A_all, f)
                # Split-half: two disjoint observed subsets of f/2 each.
                half_a = rng.binomial(A_all, f / 2)
                half_b = rng.binomial(A_all - half_a, (f / 2) / (1 - f / 2))
                for thr in THR_OBS:
                    for method in ("oracle", "observed", "split"):
                        if method == "oracle" and thr != THR_OBS[0]:
                            continue     # oracle does not depend on the obs threshold
                        if method == "split":
                            vals = [estimate(A_all[s], np.stack([half_a[s], half_b[s]]), f,
                                             method, thr) for s in range(A_all.shape[0])]
                        else:
                            vals = [estimate(A_all[s], obs_all[s], f, method, thr)
                                    for s in range(A_all.shape[0])]
                        mh = np.nanmean([v[0] for v in vals])
                        lam = np.nanmean([v[1] for v in vals])
                        n_on = int(np.mean([v[2] for v in vals]))
                        rows.append(dict(world=wname, f=f,
                                         thr=thr if method != "oracle" else np.nan,
                                         method=method, m_true=m, m_hat=mh, lam=lam,
                                         n_onsets=n_on))
        print(f"done world {wname}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv("transient_observed_validation.csv", index=False)

    mt = np.array(M_GRID)
    cross = []
    for key, g in df.groupby(["world", "f", "method", "thr"], dropna=False):
        g = g.sort_values("m_true")
        mh = g.m_hat.values
        star = crossing_point(mt, mh)
        ok = np.isfinite(mh) & (mt != 1.0)
        sign_acc = (float(np.mean(np.sign(mh[ok] - 1) == np.sign(mt[ok] - 1)))
                    if ok.any() else np.nan)
        mono = bool(np.all(np.diff(mh[np.isfinite(mh)]) > -0.02))
        cross.append(dict(world=key[0], f=key[1], method=key[2], thr=key[3], crossing=star,
                          sign_accuracy=sign_acc, monotone=mono,
                          min_onsets=int(g.n_onsets.min())))
    cr = pd.DataFrame(cross)
    cr.to_csv("transient_observed_crossings.csv", index=False)

    print("\n=== Crossing m* (ideal 1.000) and sign accuracy, by protocol ===")
    summ = cr.groupby(["method", "thr", "f"], dropna=False).agg(
        cross_min=("crossing", "min"), cross_max=("crossing", "max"),
        n_nan=("crossing", lambda s: int(s.isna().sum())),
        sign_acc_min=("sign_accuracy", "min"), all_monotone=("monotone", "all"),
        min_onsets=("min_onsets", "min"))
    print(summ.round(3).to_string())
    for meth in ("oracle", "observed", "split"):
        s = cr[cr.method == meth].crossing.dropna()
        if len(s):
            print(f"\n{meth:9}: crossing {s.min():.3f}-{s.max():.3f}, "
                  f"max |m*-1| = {np.max(np.abs(s - 1)):.3f}, "
                  f"undefined in {int(cr[cr.method == meth].crossing.isna().sum())} conditions")


if __name__ == "__main__":
    main()
