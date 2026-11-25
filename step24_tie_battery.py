"""Tie/sparsity battery for culture-arm complexity measures (A9 §9.2.1).

Implements the generators and candidates fixed in docs/A9_culture_arm_preregistration.md
§9.2.1 on 2026-11-24, before this script was written. Do not change constants without
logging an amendment.

Candidates
----------
(a) LZ:    LZ_norm = LZ76 phrase count of the integer count sequence (no binarization)
           divided by the mean phrase count of K within-series shuffles. Null value 1.
           Phrase counts use antropy.lziv_complexity(normalize=False), verified identical
           to step2's inline lz76 on multi-symbol Poisson sequences (6/6 cases, 2026-11-24).
(b) MPR_C: dC = MPR_C(x) - mean MPR_C over K shuffles; MPR_C as in the EEG arm (d = 6,
           delay 1, stable-sort tie ranking). Null value 0.
(c) MPR_C on a Gaussian-smoothed rate (sigma = 10 ms); evaluated only if (b) fails.
    A9 lists (c) for "both measures", but LZ on a continuous smoothed rate needs a
    binarization rule the spec does not give, so the LZ version of (c) is not computed
    (logged in docs/decision_log.md). If (a) fails, the LZ contrast is dropped per A9 §4.

Pass rules (A9 §9.2.1)
----------------------
LZ (a):    every iid series |LZ_norm - 1| <= 0.02; every periodic and rate-modulated series
           LZ_norm <= LZ_norm(iid at the same mean) - 0.2.
MPR_C (b): every iid series |dC| <= 0.02; every rate-modulated and bursting series
           dC > 0.02.
Series with < 2 distinct values or < 100 non-zero bins are undefined and never scored.

Usage
-----
    python3 step24_tie_battery.py
Writes tie_battery_results.csv.
"""
from __future__ import annotations

import antropy as ant
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

import step2_batch_extract as s2
from step11_param_robustness_sweep import BASE, simulate_vec

SEED = 20261124
N_BINS = 360_000                 # 1 h at 10 ms
K_SHUFFLE = 20
IID_MEANS = [0.005, 0.01, 0.02, 0.1, 0.3, 0.65, 0.75, 1.0, 1.5, 3.0, 5.0]
PERIODS = [50, 100, 200]
MOD_MEANS = [0.3, 1.0, 3.0]
MOD_DEPTH, MOD_PERIOD = 0.9, 200
BURST_F = [0.02, 0.05, 0.10]
BURST_STEPS = 3_600_000          # 1 ms steps -> 360,000 bins of 10 ms
SMOOTH_SIGMA_MS = 10.0
LZ_IID_TOL, LZ_MARGIN = 0.02, 0.2
DC_IID_TOL, DC_STRUCT = 0.02, 0.02
MIN_NONZERO = 100


def lz_count(x: np.ndarray) -> int:
    """LZ76 phrase count of an integer sequence (no binarization)."""
    return int(ant.lziv_complexity(np.asarray(x, dtype=np.int64), normalize=False))


def is_defined(x: np.ndarray) -> bool:
    """A series is scored only with >= 2 distinct values and >= MIN_NONZERO non-zero bins."""
    return len(np.unique(x)) >= 2 and int(np.count_nonzero(x)) >= MIN_NONZERO


def surrogate_stats(x: np.ndarray, rng: np.random.Generator, lz: bool = True
                    ) -> tuple[float, float]:
    """(LZ_norm, dC) of x against K within-series shuffles.

    Args:
        x: series (integer counts, or a smoothed rate for candidate (c)).
        rng: random generator for the shuffles.
        lz: compute LZ_norm (skipped for continuous smoothed input).

    Returns:
        (LZ_norm or NaN, dC)
    """
    shuf = [rng.permutation(x) for _ in range(K_SHUFFLE)]
    c_x = s2.mpr_complexity(np.asarray(x, dtype=float))[1]
    c_s = float(np.mean([s2.mpr_complexity(np.asarray(s, dtype=float))[1] for s in shuf]))
    lz_norm = np.nan
    if lz:
        lz_norm = lz_count(x) / float(np.mean([lz_count(s) for s in shuf]))
    return lz_norm, c_x - c_s


def smooth_counts(x: np.ndarray, bin_ms: float) -> np.ndarray:
    """Gaussian-smoothed rate (sigma = 10 ms) sampled at 10 ms bins."""
    sm = gaussian_filter1d(x.astype(float), SMOOTH_SIGMA_MS / bin_ms)
    return sm[:: int(round(10.0 / bin_ms))][:N_BINS]


def signals(rng: np.random.Generator) -> list[dict]:
    """All battery series: 10 ms counts, plus 1 ms counts where they exist (bursting)."""
    out = [dict(signal="silent", param=0.0, x=np.zeros(N_BINS, dtype=np.int64)),
           dict(signal="near_silent", param=0.001, x=rng.poisson(0.001, N_BINS))]
    out += [dict(signal="iid", param=mu, x=rng.poisson(mu, N_BINS)) for mu in IID_MEANS]
    for p in PERIODS:
        x = np.zeros(N_BINS, dtype=np.int64)
        x[::p] = 1
        out.append(dict(signal="periodic", param=1.0 / p, x=x))
    t = np.arange(N_BINS)
    for mu in MOD_MEANS:
        lam = mu * (1 + MOD_DEPTH * np.sin(2 * np.pi * t / MOD_PERIOD))
        out.append(dict(signal="rate_mod", param=mu, x=rng.poisson(lam)))
    A = simulate_vec(1.0, **BASE, n_seeds=1, t_steps=BURST_STEPS, seed0=SEED)[0]
    for f in BURST_F:
        fine = rng.binomial(A, f)
        x10 = fine[: N_BINS * 10].reshape(N_BINS, 10).sum(1)
        out.append(dict(signal="burst", param=f, x=x10, fine=fine))
    return out


def evaluate(df: pd.DataFrame, lz_col: str | None, dc_col: str) -> dict:
    """Apply the §9.2.1 pass rules to one candidate's columns."""
    d = df[df.defined]
    iid = d[d.signal == "iid"].set_index("param")
    out = {}
    if lz_col:
        iid_ok = bool((iid[lz_col] - 1).abs().max() <= LZ_IID_TOL)
        struct = []
        for _, r in d[d.signal.isin(["periodic", "rate_mod"])].iterrows():
            ref = iid[lz_col].iloc[int(np.abs(iid.index.values - r.param).argmin())]
            struct.append(bool(r[lz_col] <= ref - LZ_MARGIN))
        out.update(lz_pass=iid_ok and all(struct), lz_iid_ok=iid_ok,
                   lz_struct=f"{sum(struct)}/{len(struct)}")
    iid_ok = bool(iid[dc_col].abs().max() <= DC_IID_TOL)
    struct = [bool(r[dc_col] > DC_STRUCT)
              for _, r in d[d.signal.isin(["rate_mod", "burst"])].iterrows()]
    out.update(mpr_pass=iid_ok and all(struct), mpr_iid_ok=iid_ok,
               mpr_struct=f"{sum(struct)}/{len(struct)}")
    return out


def main() -> None:
    rng = np.random.default_rng(SEED)
    rows = []
    for sig in signals(rng):
        x = sig["x"]
        defined = is_defined(x)
        row = dict(signal=sig["signal"], param=sig["param"], mean=float(x.mean()),
                   n_nonzero=int(np.count_nonzero(x)), defined=defined,
                   lz_norm=np.nan, dC=np.nan, dC_smooth=np.nan)
        if defined:
            row["lz_norm"], row["dC"] = surrogate_stats(x, rng)
            sm = smooth_counts(sig["fine"], 1.0) if "fine" in sig else smooth_counts(x, 10.0)
            _, row["dC_smooth"] = surrogate_stats(sm, rng, lz=False)
        rows.append(row)
        print(f"{row['signal']:>11} {row['param']:>7.3f}  defined={defined!s:5}  "
              f"LZ_norm={row['lz_norm']:.3f}  dC={row['dC']:+.4f}  "
              f"dC_smooth={row['dC_smooth']:+.4f}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv("tie_battery_results.csv", index=False)

    a = evaluate(df, "lz_norm", "dC")
    c = evaluate(df, None, "dC_smooth")
    print("\n=== A9 §9.2.1 verdicts ===")
    print(f"LZ candidate (a):    {'PASS' if a['lz_pass'] else 'FAIL'}  (iid within "
          f"±{LZ_IID_TOL}: {a['lz_iid_ok']}; structure ≥ {LZ_MARGIN} below iid: {a['lz_struct']})")
    print(f"MPR_C candidate (b): {'PASS' if a['mpr_pass'] else 'FAIL'}  (iid |dC| ≤ "
          f"{DC_IID_TOL}: {a['mpr_iid_ok']}; structure dC > {DC_STRUCT}: {a['mpr_struct']})")
    if not a["mpr_pass"]:
        print(f"MPR_C candidate (c, smoothed): {'PASS' if c['mpr_pass'] else 'FAIL'}  "
              f"(iid: {c['mpr_iid_ok']}; structure: {c['mpr_struct']})")
    print("Silent series flagged undefined:", bool(not df[df.signal == "silent"].defined.iloc[0]))


if __name__ == "__main__":
    main()
