
from __future__ import annotations

import sys
import os
import glob
import logging
from typing import Iterator, Sequence

import numpy as np
import pandas as pd
import h5py
from scipy.ndimage import median_filter

logging.getLogger("mrestimator").setLevel(logging.ERROR)




CACHE_DIR: str = os.path.expanduser("~/allen_cache")


RUN_THR: float = 3.0
STAT_THR: float = 0.5
KMAX_BINS: int = 60
MIN_TAU_BINS: float = 2.0
MAX_TAU_FRAC: float = 0.8
M_CEIL_SCAN: float = 1.5
MIN_UNITS: int = 20
AMP_CUTOFF_MAX: float = 0.1
PRESENCE_MIN: float = 0.9
ISI_VIOL_MAX: float = 0.5
BIN_SWEEP_MS: tuple[int, ...] = (4, 8, 16)


SMOOTH_WIN_S: float = 1.0
T_TRIAL_S: float = 5.0
MAX_GAP_S: float = 0.1

SESSION_GLOB: str = "session_*.nwb"
EXPECTED_BYTES: int = 2_500_000_000

SPONT_START = "intervals/spontaneous_presentations/start_time"
SPONT_STOP = "intervals/spontaneous_presentations/stop_time"
RUN_DATA = "processing/running/running_speed/data"
RUN_TIME = "processing/running/running_speed/timestamps"
ELEC_ID = "general/extracellular_ephys/electrodes/id"
ELEC_LOC = "general/extracellular_ephys/electrodes/location"


RUNS: dict[str, tuple] = {
    "pooled":      ("pooled",     None, "none",   T_TRIAL_S, BIN_SWEEP_MS),
    "primary":     ("stationary", 1.0,  "mean",   T_TRIAL_S, BIN_SWEEP_MS),
    "smooth0.5":   ("stationary", 0.5,  "mean",   T_TRIAL_S, (4,)),
    "smooth2.0":   ("stationary", 2.0,  "mean",   T_TRIAL_S, (4,)),
    "median1.0":   ("stationary", 1.0,  "median", T_TRIAL_S, (4,)),
    "trial2s":     ("stationary", 1.0,  "mean",   2.0,       (4,)),
    "trial10s":    ("stationary", 1.0,  "mean",   10.0,      (4,)),
    "running":     ("running",    1.0,  "mean",   T_TRIAL_S, (4,)),
}


def decode(arr: np.ndarray) -> np.ndarray:
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def session_files(cache_dir: str) -> Iterator[tuple[str, str]]:
    for path in sorted(glob.glob(os.path.join(cache_dir, SESSION_GLOB))):
        sid = os.path.basename(path).replace("session_", "").replace(".nwb", "")
        if os.path.getsize(path) < EXPECTED_BYTES:
            print(f"  SKIP {sid}: truncated")
            continue
        yield sid, path


def load_units(f: h5py.File) -> pd.DataFrame:
    qual = decode(f["units/quality"][:])
    amp = np.asarray(f["units/amplitude_cutoff"][:], dtype=float)
    pres = np.asarray(f["units/presence_ratio"][:], dtype=float)
    isi = np.asarray(f["units/isi_violations"][:], dtype=float)
    peak_ch = np.asarray(f["units/peak_channel_id"][:], dtype=np.int64)

    e_id = np.asarray(f[ELEC_ID][:], dtype=np.int64)
    e_loc = decode(f[ELEC_LOC][:])
    loc_map = dict(zip(e_id.tolist(), e_loc.tolist()))
    region = np.array([loc_map.get(int(c), "unknown") for c in peak_ch])

    good = (qual == "good") & (amp <= AMP_CUTOFF_MAX) & \
           (pres >= PRESENCE_MIN) & (isi <= ISI_VIOL_MAX)
    return pd.DataFrame({"row": np.arange(len(qual)), "region": region, "good": good})


def unit_spikes(f: h5py.File) -> tuple[np.ndarray, np.ndarray]:
    st = np.asarray(f["units/spike_times"][:], dtype=float)
    idx = np.asarray(f["units/spike_times_index"][:], dtype=np.int64)
    starts = np.r_[0, idx[:-1]]
    return st, np.c_[starts, idx]


def smooth_signed(speed: np.ndarray, win_n: int, kind: str) -> np.ndarray:
    if win_n <= 1 or kind == "none":
        return speed
    if kind == "median":
        return median_filter(speed, size=win_n, mode="nearest")
    return np.convolve(speed, np.ones(win_n) / win_n, mode="same")


def in_spont(ts: np.ndarray, starts: np.ndarray, stops: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    o = np.argsort(starts)
    s0, s1 = starts[o], stops[o]
    assert np.all(s0[1:] >= s1[:-1] - 1e-9), "spontaneous intervals overlap"
    idx = np.searchsorted(s0, ts, side="right") - 1
    valid = idx >= 0
    mask = np.zeros(ts.shape, dtype=bool)
    mask[valid] = ts[valid] < s1[idx[valid]]
    return mask, idx


def state_bouts(ts: np.ndarray, keep: np.ndarray, iid: np.ndarray) -> list[tuple[float, float]]:
    if len(ts) == 0:
        return []
    flip = np.r_[True, (keep[1:] != keep[:-1]) | (iid[1:] != iid[:-1])]
    grp = np.cumsum(flip) - 1
    out = []
    for g in np.unique(grp[keep]):
        sel = np.where(grp == g)[0]
        out.append((float(ts[sel[0]]), float(ts[sel[-1]])))
    return out


def build_windows(path: str, partition: str, win_s: float | None,
                  kind: str) -> tuple[list[tuple[float, float]], dict]:
    with h5py.File(path, "r") as f:
        starts = np.asarray(f[SPONT_START][:], dtype=float)
        stops = np.asarray(f[SPONT_STOP][:], dtype=float)
        if partition == "pooled":
            return list(zip(starts.tolist(), stops.tolist())), {}
        speed = np.asarray(f[RUN_DATA][:], dtype=float)
        ts = np.asarray(f[RUN_TIME][:], dtype=float)

    o = np.argsort(ts)
    ts, speed = ts[o], speed[o]
    ok = np.isfinite(ts) & np.isfinite(speed)
    ts, speed = ts[ok], speed[ok]

    fs = 1.0 / float(np.median(np.diff(ts)))
    win_n = max(1, int(round((win_s or 0.0) * fs)))
    sm = smooth_signed(speed, win_n, kind)
    mag = np.abs(sm)

    mask, idx_all = in_spont(ts, starts, stops)
    ts_s, mag_s, iid = ts[mask], mag[mask], idx_all[mask]

    keep = (mag_s <= STAT_THR) if partition == "stationary" else (mag_s >= RUN_THR)
    bouts = state_bouts(ts_s, keep, iid)

    resid = np.abs(speed - sm)
    diag = {"fs_hz": round(fs, 2), "resid_sd_all": float(resid.std())}
    return bouts, diag


def segment(windows: Sequence[tuple[float, float]], t_trial: float) -> list[tuple[float, float]]:
    trials = []
    for t0, t1 in windows:
        n = int((t1 - t0) // t_trial)
        trials.extend((t0 + k * t_trial, t0 + (k + 1) * t_trial) for k in range(n))
    return trials


def estimate_m(counts: np.ndarray) -> dict:
    import mrestimator as mre
    out = dict(m_MR=np.nan, tau_bins=np.nan, mr_ok=False, reason="")
    if counts.shape[0] < 3:
        out["reason"] = "too few trials"; return out
    if counts.sum() < 20:
        out["reason"] = "too few events"; return out
    try:
        rk = mre.coefficients(counts, steps=(1, KMAX_BINS), dt=1, method="trialseparated")
        ft = mre.fit(rk, fitfunc="exponential")
        m, tau = float(ft.mre), float(ft.tau)
    except Exception as e:
        out["reason"] = f"fit failed: {type(e).__name__}"; return out
    out.update(m_MR=m, tau_bins=tau)
    if not np.isfinite(m) or not np.isfinite(tau):
        out["reason"] = "non-finite"
    elif tau < MIN_TAU_BINS:
        out["reason"] = "degenerate (tau too small)"
    elif tau > MAX_TAU_FRAC * KMAX_BINS:
        out["reason"] = "no decay within window"
    elif not (0 < m < M_CEIL_SCAN):
        out["reason"] = "m out of physical range"
    else:
        out["mr_ok"] = True
    return out


def scan_session(sid: str, path: str, label: str) -> list[dict]:
    partition, win_s, kind, t_trial, bins_ms = RUNS[label]
    windows, _ = build_windows(path, partition, win_s, kind)
    trials = segment(windows, t_trial)
    if len(trials) < 3:
        print(f"    {sid}: only {len(trials)} trials — skipping")
        return []

    with h5py.File(path, "r") as f:
        units = load_units(f)
        st, offs = unit_spikes(f)
        good = units[units["good"]]
        rows = []
        for region, grp in good.groupby("region"):
            if region in ("unknown", "root", "") or len(grp) < MIN_UNITS:
                continue
            spikes = np.sort(np.concatenate(
                [st[offs[r, 0]:offs[r, 1]] for r in grp["row"].to_numpy()]))
            for bin_ms in bins_ms:
                bin_s = bin_ms / 1000.0
                n_bins = int(round(t_trial / bin_s))
                mat = np.empty((len(trials), n_bins), dtype=np.int64)
                for i, (t0, _) in enumerate(trials):
                    edges = t0 + np.arange(n_bins + 1) * bin_s
                    lo = np.searchsorted(spikes, edges[0], "left")
                    hi = np.searchsorted(spikes, edges[-1], "left")
                    mat[i] = np.histogram(spikes[lo:hi], bins=edges)[0]
                if n_bins <= KMAX_BINS + 2:
                    rows.append(dict(run=label, session_id=sid, region=region,
                                     bin_ms=bin_ms, m_MR=np.nan, mr_ok=False,
                                     reason="trial shorter than KMAX+2 bins"))
                    continue
                res = estimate_m(mat)
                rows.append(dict(
                    run=label, session_id=sid, region=region, partition=partition,
                    smooth_win_s=win_s, smooth_kind=kind, t_trial_s=t_trial,
                    bin_ms=bin_ms, n_units=int(len(grp)), n_trials=int(len(trials)),
                    total_bins=int(mat.size), mean_rate_hz=float(mat.sum() / (len(trials) * t_trial)),
                    **res))
    return rows


def residcheck(paths: list[tuple[str, str]]) -> None:
    rows = []
    for sid, path in paths:
        with h5py.File(path, "r") as f:
            starts = np.asarray(f[SPONT_START][:], dtype=float)
            stops = np.asarray(f[SPONT_STOP][:], dtype=float)
            speed = np.asarray(f[RUN_DATA][:], dtype=float)
            ts = np.asarray(f[RUN_TIME][:], dtype=float)
        o = np.argsort(ts); ts, speed = ts[o], speed[o]
        ok = np.isfinite(ts) & np.isfinite(speed); ts, speed = ts[ok], speed[ok]
        fs = 1.0 / float(np.median(np.diff(ts)))
        sm = smooth_signed(speed, max(1, int(round(SMOOTH_WIN_S * fs))), "mean")
        mask, idx_all = in_spont(ts, starts, stops)
        stat = np.abs(sm)[mask] <= STAT_THR
        r_all = np.abs(speed - sm)
        r_stat = r_all[mask][stat]
        rows.append({
            "session_id": sid, "fs_hz": round(fs, 2),
            "resid_sd_all": float(r_all.std()),
            "resid_p99_all": float(np.percentile(r_all, 99)),
            "resid_sd_stationary": float(r_stat.std()),
            "resid_p99_stationary": float(np.percentile(r_stat, 99)),
            "glitch_flag_stationary": bool(np.percentile(r_stat, 99) > 3 * r_stat.std()),
        })
    df = pd.DataFrame(rows)
    df.to_csv("allen_m_scan_stationary_residcheck.csv", index=False)
    print("\nRESIDUAL CHECK — restricted to stationary bouts (win=1.0s, mean)")
    print(df.to_string(index=False))
    print("\n  If resid_sd_stationary << resid_sd_all, step14b's glitch criterion fired")
    print("  on locomotion transients, not glitches, and the mean-smoothing choice stands.")
    print("  The median1.0 run is the independent check either way.")


def main() -> None:
    args = sys.argv[1:]
    cache = next((a for a in args if not a.startswith("--")), CACHE_DIR)
    cache = os.path.expanduser(cache)
    want = next((a.split("=", 1)[-1] for a in args if a.startswith("--runs")), None)
    if want and "=" not in want:
        want = args[args.index("--runs") + 1] if "--runs" in args else None
    labels = want.split(",") if want else list(RUNS)

    paths = list(session_files(cache))
    if not paths:
        print(f"No sessions in {cache}"); return

    if "--residcheck" in args:
        residcheck(paths)
        if not want:
            labels = list(RUNS)

    all_rows = []
    for label in labels:
        if label not in RUNS:
            print(f"  unknown run '{label}' — skipping"); continue
        print(f"\n=== RUN: {label}  {RUNS[label][:4]} ===")
        for sid, path in paths:
            try:
                rows = scan_session(sid, path, label)
            except (KeyError, ValueError, OSError, AssertionError) as e:
                print(f"    FAIL {sid}: {type(e).__name__}: {e}"); continue
            all_rows.extend(rows)
            ok = [r for r in rows if r.get("mr_ok")]
            if ok:
                mx = max(r["m_MR"] for r in ok)
                print(f"    {sid}: {len(ok):3d}/{len(rows):3d} mr_ok | max m_MR {mx:.3f} | "
                      f"n>=0.98: {sum(r['m_MR'] >= 0.98 for r in ok)} | "
                      f">1: {sum(r['m_MR'] > 1.0 for r in ok)}")
            else:
                print(f"    {sid}: 0/{len(rows)} mr_ok")

    if not all_rows:
        print("\nNothing scanned."); return
    df = pd.DataFrame(all_rows)
    df.to_csv("allen_m_scan_stationary.csv", index=False)

    print("\n" + "=" * 78)
    print("SUMMARY BY RUN x BIN  (criterion 5 gate = 15% mr_ok)")
    print("=" * 78)
    for (run, bin_ms), g in df.groupby(["run", "bin_ms"]):
        ok = g[g["mr_ok"] == True]
        frac = len(ok) / len(g) if len(g) else 0.0
        if len(ok):
            print(f"  {run:11s} bin {bin_ms:2d}ms | cells {len(g):3d} | mr_ok {frac:5.1%} | "
                  f"max m {ok['m_MR'].max():.3f} | median m {ok['m_MR'].median():.3f} | "
                  f">=0.98: {int((ok['m_MR'] >= 0.98).sum())} | >1: {int((ok['m_MR'] > 1.0).sum())}")
        else:
            print(f"  {run:11s} bin {bin_ms:2d}ms | cells {len(g):3d} | mr_ok {frac:5.1%} | —")

    if "pooled" in df["run"].values:
        print("\n" + "-" * 78)
        print("CONCORDANCE GATE — compare `pooled` against allen_m_scan_binsweep.csv:")
        print("  known: max m_MR <= 0.977 at every bin width; 0 of 33-44 cells >= 0.98;")
        print("         0 supercritical.")
        print("  Concordant -> estimator path verified; stationary rows are interpretable.")
        print("  Discordant -> STOP. Reconcile against the real step7 source first.")
        print("                Do not explain a discordance. Every stationary number")
        print("                is uninterpretable until the paths agree.")
        print("-" * 78)

    print("\nmr_ok drop reasons (criterion 5 reporting):")
    bad = df[df["mr_ok"] == False]
    if len(bad):
        print(bad.groupby(["run", "reason"]).size().to_string())
    else:
        print("  none")


if __name__ == "__main__":
    main()
