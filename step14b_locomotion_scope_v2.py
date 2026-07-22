
from __future__ import annotations

import sys
import glob
import os
from typing import Iterator

import numpy as np
import pandas as pd
import h5py




CACHE_DIR: str = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/allen_cache")

RUN_THR: float = 3.0
STAT_THR: float = 0.5
MIN_BOUT_S: float = 0.25
N_BOUTS_MIN: int = 3
MIN_POOLED_STAT_S: float = 60.0
MAX_GAP_S: float = 0.1


SMOOTH_WINDOWS_S: tuple[float, ...] = (0.0, 0.25, 0.5, 1.0, 2.0)
REPORT_WINDOW_S: float = 1.0

SESSION_GLOB: str = "session_*.nwb"
EXPECTED_BYTES: int = 2_500_000_000

SPONT_START = "intervals/spontaneous_presentations/start_time"
SPONT_STOP = "intervals/spontaneous_presentations/stop_time"
RUN_DATA = "processing/running/running_speed/data"
RUN_TIME = "processing/running/running_speed/timestamps"


def session_files(cache_dir: str) -> Iterator[tuple[str, str]]:
    for path in sorted(glob.glob(os.path.join(cache_dir, SESSION_GLOB))):
        sid = os.path.basename(path).replace("session_", "").replace(".nwb", "")
        size = os.path.getsize(path)
        if size < EXPECTED_BYTES:
            print(f"  SKIP {sid}: {size/1e9:.2f} GB — looks truncated")
            continue
        yield sid, path


def load_running(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as f:
        for p in (SPONT_START, SPONT_STOP, RUN_DATA, RUN_TIME):
            if p not in f:
                raise KeyError(f"missing HDF5 path {p}")
        starts = np.asarray(f[SPONT_START][:], dtype=float)
        stops = np.asarray(f[SPONT_STOP][:], dtype=float)
        speed = np.asarray(f[RUN_DATA][:], dtype=float)
        ts = np.asarray(f[RUN_TIME][:], dtype=float)
    if len(speed) != len(ts):
        raise ValueError(f"running data/timestamps mismatch ({len(speed)} vs {len(ts)})")
    order = np.argsort(ts)
    ts, speed = ts[order], speed[order]
    ok = np.isfinite(ts) & np.isfinite(speed)
    return ts[ok], speed[ok], starts, stops


def smooth_signed(speed: np.ndarray, win_samples: int) -> np.ndarray:
    if win_samples <= 1:
        return speed
    kern = np.ones(win_samples, dtype=float) / win_samples
    return np.convolve(speed, kern, mode="same")


def in_spontaneous(ts: np.ndarray, starts: np.ndarray, stops: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(starts)
    s0, s1 = starts[order], stops[order]
    assert np.all(s0[1:] >= s1[:-1] - 1e-9), "spontaneous intervals overlap"
    idx = np.searchsorted(s0, ts, side="right") - 1
    valid = idx >= 0
    mask = np.zeros(ts.shape, dtype=bool)
    mask[valid] = ts[valid] < s1[idx[valid]]
    return mask, idx


def stationary_bouts(ts: np.ndarray, stat: np.ndarray, iid: np.ndarray) -> pd.DataFrame:
    if len(ts) == 0:
        return pd.DataFrame(columns=["interval_id", "t_start", "t_stop", "duration_s", "n_samples"])
    flip = np.r_[True, (stat[1:] != stat[:-1]) | (iid[1:] != iid[:-1])]
    grp = np.cumsum(flip) - 1
    rows = []
    for g in np.unique(grp[stat]):
        sel = np.where(grp == g)[0]
        rows.append({
            "interval_id": int(iid[sel[0]]),
            "t_start": float(ts[sel[0]]),
            "t_stop": float(ts[sel[-1]]),
            "duration_s": float(ts[sel[-1]] - ts[sel[0]]),
            "n_samples": int(len(sel)),
        })
    return pd.DataFrame(rows, columns=["interval_id", "t_start", "t_stop",
                                       "duration_s", "n_samples"])


def scan(sid: str, path: str) -> tuple[list[dict], pd.DataFrame]:
    ts, speed, starts, stops = load_running(path)
    dt_med = float(np.median(np.diff(ts)))
    fs = 1.0 / dt_med

    rows, report_bouts = [], pd.DataFrame()
    for win_s in SMOOTH_WINDOWS_S:
        win_n = max(1, int(round(win_s * fs)))


        sm = smooth_signed(speed, win_n)
        mag = np.abs(sm)

        mask, idx_all = in_spontaneous(ts, starts, stops)
        ts_s, mag_s, iid = ts[mask], mag[mask], idx_all[mask]
        if len(ts_s) < 2:
            raise ValueError("no running samples inside spontaneous intervals")

        w = np.clip(np.gradient(ts_s), 0.0, MAX_GAP_S)
        total = float(w.sum())
        running = mag_s >= RUN_THR
        stat = mag_s <= STAT_THR
        mid = ~running & ~stat

        bouts = stationary_bouts(ts_s, stat, iid)
        fittable = bouts[bouts["duration_s"] >= MIN_BOUT_S]

        resid = np.abs(speed - sm)
        rows.append({
            "session_id": sid,
            "smooth_win_s": win_s,
            "fs_hz": round(fs, 2),
            "spont_covered_s": total,
            "frac_running": float(w[running].sum() / total),
            "frac_stationary": float(w[stat].sum() / total),
            "frac_intermediate": float(w[mid].sum() / total),
            "stationary_s": float(w[stat].sum()),
            "n_stat_bouts": int(len(bouts)),
            "n_stat_bouts_fittable": int(len(fittable)),
            "stationary_fittable_s": float(fittable["duration_s"].sum()),
            "longest_stat_bout_s": float(bouts["duration_s"].max()) if len(bouts) else 0.0,
            "median_speed_cms": float(np.median(mag_s)),
            "resid_sd": float(resid.std()),
            "resid_p99": float(np.percentile(resid, 99)),
        })
        if win_s == REPORT_WINDOW_S:
            report_bouts = bouts.assign(session_id=sid)
    return rows, report_bouts


def main() -> None:
    print(f"Locomotion scope v2 (rectification bug fixed) — cache: {CACHE_DIR}")
    print(f"RUN >= {RUN_THR} | STAT <= {STAT_THR} cm/s (UNCHANGED) | "
          f"smoothing sweep: {SMOOTH_WINDOWS_S}\n")

    all_rows, all_bouts = [], []
    for sid, path in session_files(CACHE_DIR):
        try:
            rows, bouts = scan(sid, path)
        except (KeyError, ValueError, OSError) as e:
            print(f"  FAIL {sid}: {type(e).__name__}: {e}")
            continue
        all_rows.extend(rows)
        if len(bouts):
            all_bouts.append(bouts)
        print(f"  {sid}:")
        for r in rows:
            tag = " (v1 reference)" if r["smooth_win_s"] == 0.0 else ""
            print(f"    win {r['smooth_win_s']:4.2f}s | run {r['frac_running']:5.1%} | "
                  f"stat {r['frac_stationary']:5.1%} | mid {r['frac_intermediate']:5.1%} | "
                  f"fittable bouts {r['n_stat_bouts_fittable']:4d} "
                  f"({r['stationary_fittable_s']:7.1f}s) | "
                  f"longest {r['longest_stat_bout_s']:6.2f}s{tag}")

    if not all_rows:
        print("\nNo sessions scanned.")
        return

    df = pd.DataFrame(all_rows)
    df.to_csv("locomotion_scope_v2_sweep.csv", index=False)
    if all_bouts:
        pd.concat(all_bouts, ignore_index=True).to_csv("locomotion_scope_v2_bouts.csv", index=False)

    print("\n" + "=" * 78)
    print("POOLED BY SMOOTHING WINDOW (time-weighted)")
    print("=" * 78)
    for win_s in SMOOTH_WINDOWS_S:
        sub = df[df["smooth_win_s"] == win_s]
        tw = sub["spont_covered_s"]
        p_run = float((sub["frac_running"] * tw).sum() / tw.sum())
        p_stat = float((sub["frac_stationary"] * tw).sum() / tw.sum())
        p_mid = float((sub["frac_intermediate"] * tw).sum() / tw.sum())
        n_fit = int(sub["n_stat_bouts_fittable"].sum())
        s_fit = float(sub["stationary_fittable_s"].sum())
        viable = (n_fit >= N_BOUTS_MIN) and (s_fit >= MIN_POOLED_STAT_S)
        tag = " (v1 reference — known biased)" if win_s == 0.0 else ""
        print(f"  win {win_s:4.2f}s | run {p_run:5.1%} | stat {p_stat:5.1%} | "
              f"mid {p_mid:5.1%} | fittable {n_fit:5d} bouts / {s_fit:8.1f}s | "
              f"{'VIABLE' if viable else 'thin  '}{tag}")

    print("\nSTABILITY CHECK — the verdict must hold across windows >= 0.25s.")
    print("If it flips window to window, that instability IS the result. Report it;")
    print("do not select the window that yields the most stationary data.\n")
    print("Noise diagnostic (raw minus smoothed, at win=1.0s):")
    print(df[df["smooth_win_s"] == 1.0][
        ["session_id", "fs_hz", "resid_sd", "resid_p99", "median_speed_cms"]
    ].to_string(index=False))
    print("\n  resid_sd is the dither amplitude that broke v1. If resid_p99 >> 3*resid_sd,")
    print("  the trace carries isolated glitches and a median-smoothing variant is the")
    print("  logged fallback (a separate decision, not a silent swap).")


if __name__ == "__main__":
    main()
