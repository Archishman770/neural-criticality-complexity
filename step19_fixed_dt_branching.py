"""Branching-ratio estimation with a FIXED bin width, as a confound control for step3b.

Motivation
----------
`step3b_branching_ratio_v2.py` sets the bin width per cell from the mean
inter-event interval (`choose_dt`, dt = mean IEI). Because dt is therefore
1/event-rate, and because both Lempel-Ziv and MPR_C rise with event rate, the
reported complexity-criticality relationship (manuscript s3.1) has dt on both
sides of it. This script re-estimates the branching ratio holding dt constant
across all subject-condition cells, so that any surviving relationship cannot be
carried by the binning parameter.

Everything other than dt selection is identical to step3b: same 2.5-SD peak
detection, same trial-separated MR coefficients, same 60-bin fit window, same
quality criteria.

Usage
-----
    python3 step19_fixed_dt_branching.py <data_dir>

Writes `branching_summary_fixeddt.csv` (long format: one row per cell per fixed
dt) into <data_dir>.
"""

from __future__ import annotations

import glob
import logging
import os
import sys
import time

import numpy as np
from scipy.signal import find_peaks

logging.getLogger("mrestimator").setLevel(logging.ERROR)

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
THRESHOLD = 2.5
KMAX_BINS = 60
MIN_TAU_BINS = 2.0
MAX_TAU_FRAC = 0.8
NEG_ONLY = False
LEVELS = {1: "baseline", 2: "mild", 3: "moderate", 4: "recovery"}

# The dt values observed in branching_summary_v2.csv span 5-11 samples.
# Sweeping the whole range shows whether conclusions depend on which fixed
# value is chosen.
FIXED_DTS = (5, 6, 7, 8, 9, 10, 11)


def norm(name: str) -> str:
    """Normalise a filename stem for matching against datainfo.mat keys."""
    return name.replace(" ", "_").replace(".", "_")


def population_activity(
    epoch: np.ndarray, thr: float = THRESHOLD, neg_only: bool = NEG_ONLY
) -> np.ndarray:
    """Per-sample count of supra-threshold events across channels.

    Args:
        epoch: (n_channels, n_samples) array for a single epoch.
        thr: peak height threshold in units of per-channel SD.
        neg_only: if True use negative deflections only, else absolute value.

    Returns:
        (n_samples,) array of event counts.
    """
    z = (epoch - epoch.mean(1, keepdims=True)) / (epoch.std(1, keepdims=True) + 1e-12)
    a = np.zeros(z.shape[1])
    for ch in range(z.shape[0]):
        sig = -z[ch] if neg_only else np.abs(z[ch])
        peaks, _ = find_peaks(sig, height=thr)
        a[peaks] += 1
    return a


def rebin(a: np.ndarray, dt: int) -> np.ndarray:
    """Sum an activity trace into non-overlapping bins of width dt samples."""
    nb = len(a) // dt
    return a[: nb * dt].reshape(nb, dt).sum(1) if nb else a[:0]


def estimate_m_at_dt(activities: list[np.ndarray], fs: float, dt: int) -> dict:
    """Multistep-regression branching ratio at a caller-specified bin width.

    Identical to step3b's `estimate_m` except that dt is supplied rather than
    derived from the data.

    Args:
        activities: list of per-epoch population activity traces.
        fs: sampling frequency in Hz.
        dt: bin width in samples.

    Returns:
        dict with dt, m_MR, tau_bins, tau_ms, mr_ok, reason.
    """
    import mrestimator as mre

    binned = [rebin(a, dt) for a in activities]
    binned = [b for b in binned if len(b) > KMAX_BINS + 2]
    out = dict(dt=dt, m_MR=np.nan, tau_bins=np.nan, tau_ms=np.nan, mr_ok=False, reason="")
    if len(binned) < 3:
        out["reason"] = "too few usable epochs"
        return out
    L = min(len(b) for b in binned)
    arr = np.vstack([b[:L] for b in binned])
    if arr.sum() < 20:
        out["reason"] = "too few events"
        return out
    try:
        rk = mre.coefficients(arr, steps=(1, KMAX_BINS), dt=1, method="trialseparated")
        ft = mre.fit(rk, fitfunc="exponential")
        m, tau = float(ft.mre), float(ft.tau)
    except Exception as e:  # noqa: BLE001 - mirror step3b's behaviour
        out["reason"] = f"fit failed: {type(e).__name__}"
        return out

    out.update(m_MR=m, tau_bins=tau, tau_ms=tau * dt / fs * 1000)

    if not np.isfinite(m) or not np.isfinite(tau):
        out["reason"] = "non-finite"
    elif tau < MIN_TAU_BINS:
        out["reason"] = "degenerate (tau too small)"
    elif tau > MAX_TAU_FRAC * KMAX_BINS:
        out["reason"] = "no decay within window"
    elif not (0 < m < 1.2):
        out["reason"] = "m out of range"
    else:
        out["mr_ok"] = True
        out["reason"] = "ok"
    return out


def condition_map(datainfo_path: str) -> dict:
    """Map normalised dataset stem -> sedation level, condition, concentration."""
    from pymatreader import read_mat

    rows = read_mat(datainfo_path)["datainfo"]
    return {
        norm(r[0]): dict(level=int(r[1]), condition=LEVELS[int(r[1])], conc=r[2])
        for r in rows
    }


def main() -> None:
    import mne
    import pandas as pd

    datainfo = os.path.join(DATA_DIR, "datainfo.mat")
    if not os.path.exists(datainfo):
        sys.exit(f"datainfo.mat not found in {DATA_DIR!r}")
    cmap = condition_map(datainfo)
    set_files = sorted(glob.glob(os.path.join(DATA_DIR, "**", "*.set"), recursive=True))
    if not set_files:
        sys.exit(f"no .set files under {DATA_DIR!r}")
    print(f"Found {len(set_files)} .set files. Fixed dt sweep: {FIXED_DTS}\n")

    records = []
    t0 = time.time()
    for fi, f in enumerate(set_files, 1):
        stem = norm(os.path.basename(f)[:-4])
        meta = cmap.get(stem)
        if meta is None:
            print(f"  [{fi}/{len(set_files)}] {stem}: not in datainfo, skipping")
            continue
        subject = os.path.basename(f).split("-")[0]
        ep = mne.io.read_epochs_eeglab(f, verbose="ERROR")
        fs = ep.info["sfreq"]
        X = ep.get_data()
        # Peak detection is the expensive step; do it once and reuse per dt.
        acts = [population_activity(X[e]) for e in range(len(X))]
        rate = float(np.mean([a.sum() for a in acts]))

        msgs = []
        for dt in FIXED_DTS:
            r = estimate_m_at_dt(acts, fs, dt)
            records.append(
                dict(
                    subject=subject,
                    level=meta["level"],
                    condition=meta["condition"],
                    conc=meta["conc"],
                    dt_fixed=dt,
                    m_MR=r["m_MR"],
                    tau_ms=r["tau_ms"],
                    tau_bins=r["tau_bins"],
                    mr_ok=r["mr_ok"],
                    reason=r["reason"],
                    events_per_epoch=rate,
                    n_epochs=len(X),
                )
            )
            msgs.append(f"{dt}:{r['m_MR']:.3f}{'' if r['mr_ok'] else '*'}")
        print(
            f"  [{fi}/{len(set_files)}] {stem[:26]:26} {meta['condition']:8} "
            f"{' '.join(msgs)} ({time.time()-t0:.0f}s)"
        )
        del X, ep, acts

    df = pd.DataFrame.from_records(records)
    out = os.path.join(DATA_DIR, "branching_summary_fixeddt.csv")
    df.to_csv(out, index=False)
    print(f"\nDone in {time.time()-t0:.0f}s. {len(df)} rows -> {out}")
    print(df.groupby("dt_fixed").agg(n_ok=("mr_ok", "sum"), m_med=("m_MR", "median")))


if __name__ == "__main__":
    main()
