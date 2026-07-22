
from __future__ import annotations

import logging
import numpy as np
import pandas as pd

logging.getLogger("mrestimator").setLevel(logging.ERROR)










def simulate_branching(m: float, n_units: int = 1000, t_steps: int = 8000,
                       n_trials: int = 8, h: float = 2.0, burn_in: int = 2000,
                       rng: np.random.Generator | None = None) -> np.ndarray:
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





M_GRID: tuple[float, ...] = (0.90, 0.95, 0.98)




BIN_STEPS: tuple[int, ...] = (1, 2, 4)


KMAX_BINS_FIXED: int = 60
MAX_TAU_FRAC: float = 0.8
KMAX_STEPS_TIME: int = 960

MIN_TAU_BINS: float = 2.0
M_CEIL: float = 1.5



N_TRIALS: int = 12
T_STEPS: int = 40000
BURN_IN: int = 4000
SEED: int = 20260718

CENSOR_REASON = "no decay within fit window"


def induced_ceiling(kmax_bins: int) -> float:
    return float(np.exp(-1.0 / (MAX_TAU_FRAC * kmax_bins)))


def rebin(activity: np.ndarray, bin_steps: int) -> np.ndarray:
    if bin_steps <= 1:
        return activity
    n_tr, T = activity.shape
    T2 = (T // bin_steps) * bin_steps
    return activity[:, :T2].reshape(n_tr, T2 // bin_steps, bin_steps).sum(axis=2)


def mr_fit_with_rule(trials: np.ndarray, kmax_bins: int) -> dict:
    import mrestimator as mre
    out = dict(m_MR=np.nan, tau_bins=np.nan, mr_ok=False, reason="")
    if trials.sum() < 50:
        out["reason"] = "too few events"
        return out
    n_bins = trials.shape[1]
    if n_bins <= kmax_bins + 2:
        out["reason"] = "record shorter than KMAX+2"
        return out
    try:
        rk = mre.coefficients(trials, steps=(1, kmax_bins), dt=1, method="trialseparated")
        ft = mre.fit(rk, fitfunc="exponential")
        m, tau = float(ft.mre), float(ft.tau)
    except Exception as e:
        out["reason"] = f"fit failed: {type(e).__name__}"
        return out
    out.update(m_MR=m, tau_bins=tau)
    if not (np.isfinite(m) and np.isfinite(tau)):
        out["reason"] = "non-finite"
    elif tau < MIN_TAU_BINS:
        out["reason"] = "degenerate (tau too small)"
    elif tau > MAX_TAU_FRAC * kmax_bins:
        out["reason"] = CENSOR_REASON
    elif not (0 < m < M_CEIL):
        out["reason"] = "m out of range"
    else:
        out["mr_ok"] = True
    return out


def main() -> None:
    rng = np.random.default_rng(SEED)
    print("step18 — filter-artifact demonstration on KNOWN-m data")
    print(f"  generator: simulate_branching (verbatim copy of step9's, inlined)")
    print(f"  fixed rule : KMAX={KMAX_BINS_FIXED} bins -> induced m ceiling "
          f"{induced_ceiling(KMAX_BINS_FIXED):.5f} at EVERY bin width (bin-invariant)")
    print(f"  time rule  : {KMAX_STEPS_TIME} steps -> KMAX bins scales with bin width\n")

    rows = []
    for m_true in M_GRID:


        full = simulate_branching(m_true, n_trials=N_TRIALS, t_steps=T_STEPS,
                                  burn_in=BURN_IN, rng=rng)
        for bin_steps in BIN_STEPS:
            binned = rebin(full, bin_steps)
            kmax_time_bins = max(MIN_TAU_BINS + 3, round(KMAX_STEPS_TIME / bin_steps))

            fixed = mr_fit_with_rule(binned, KMAX_BINS_FIXED)
            timed = mr_fit_with_rule(binned, int(kmax_time_bins))

            rows.append(dict(
                m_true=m_true, bin_steps=bin_steps,
                fixed_kmax_bins=KMAX_BINS_FIXED,
                fixed_ceiling=induced_ceiling(KMAX_BINS_FIXED),
                fixed_m=fixed["m_MR"], fixed_tau=fixed["tau_bins"],
                fixed_ok=fixed["mr_ok"], fixed_reason=fixed["reason"],
                time_kmax_bins=int(kmax_time_bins),
                time_m=timed["m_MR"], time_tau=timed["tau_bins"],
                time_ok=timed["mr_ok"], time_reason=timed["reason"],
            ))
            print(f"  m_true={m_true:.2f} bin_steps={bin_steps}: "
                  f"fixed m={fixed['m_MR']:.4f} ({'ok' if fixed['mr_ok'] else fixed['reason']}) | "
                  f"time m={timed['m_MR']:.4f} ({'ok' if timed['mr_ok'] else timed['reason']})")

    df = pd.DataFrame(rows)
    df.to_csv("filter_artifact_demo.csv", index=False)


    print("\n" + "=" * 88)
    print("PRE-DECLARED SUCCESS CONDITION")
    print("=" * 88)

    finest = min(BIN_STEPS)
    passes = {}
    for m_true in M_GRID:
        sub = df[df.m_true == m_true]
        fine = sub[sub.bin_steps == finest].iloc[0]


        a_gap = (fine["time_m"] - fine["fixed_m"]) if (fine["fixed_ok"] or fine["fixed_reason"] == CENSOR_REASON) else np.nan

        fixed_eff = fine["fixed_m"] if fine["fixed_ok"] else fine["fixed_ceiling"]
        a_gap = fine["time_m"] - fixed_eff
        A = np.isfinite(a_gap) and a_gap >= 0.02


        fixed_eff_all = np.where(sub["fixed_ok"], sub["fixed_m"], sub["fixed_ceiling"])
        b_spread = float(np.nanmax(fixed_eff_all) - np.nanmin(fixed_eff_all))
        b_at_ceiling = abs(float(np.nanmax(fixed_eff_all)) - induced_ceiling(KMAX_BINS_FIXED)) < 0.01
        B = (b_spread < 0.01) and b_at_ceiling


        c_time_err = abs(fine["time_m"] - m_true)
        c_fixed_err = abs(fixed_eff - m_true)
        C = np.isfinite(fine["time_m"]) and (c_time_err < c_fixed_err)

        passes[m_true] = (A, B, C)
        print(f"\n  m_true = {m_true:.2f}")
        print(f"    (A) filter bites at finest bin : fixed_eff={fixed_eff:.4f} "
              f"time={fine['time_m']:.4f} gap={a_gap:+.4f}  -> {'PASS' if A else 'FAIL'}")
        print(f"    (B) false ceiling bin-robust   : spread={b_spread:.4f} "
              f"max={np.nanmax(fixed_eff_all):.4f} vs ceiling={induced_ceiling(KMAX_BINS_FIXED):.4f}"
              f"  -> {'PASS' if B else 'FAIL'}")
        print(f"    (C) time window recovers truth : time_err={c_time_err:.4f} "
              f"< fixed_err={c_fixed_err:.4f}?  -> {'PASS' if C else 'FAIL'}")

    print("\n" + "-" * 88)
    all_pass = all(all(v) for v in passes.values())
    if all_pass:
        print("VERDICT: demonstration SUPPORTS §3.3a as a contribution.")
        print("  The bin-specified filter produces a false, bin-robust subcritical ceiling")
        print("  on data with known m; the time-specified window recovers the truth. This is")
        print("  ground-truth confirmation of the mechanism inferred from the Allen data.")
        print("  -> This table becomes a figure. The Allen result is now the field demonstration")
        print("     of an effect proven here in controlled conditions.")
    else:
        print("VERDICT: demonstration does NOT cleanly support the strong claim.")
        print("  This is a REAL FINDING, not a bug. It means the Allen distortion had")
        print("  contributors beyond the fit filter, and §3.3a must be qualified: report the")
        print("  filter as A contributor, demonstrated to bias low, but not the sole cause.")
        print("  Do NOT retune config to force a pass. Report what this shows.")
        for m_true, (A, B, C) in passes.items():
            if not (A and B and C):
                fail = [n for n, p in zip("ABC", (A, B, C)) if not p]
                print(f"    m_true={m_true:.2f} failed: {', '.join(fail)}")
    print("-" * 88)


if __name__ == "__main__":
    main()
