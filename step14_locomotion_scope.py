
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
USE_ABS_SPEED: bool = True
MIN_BOUT_S: float = 0.25
N_BOUTS_MIN: int = 3
MAX_GAP_S: float = 0.1

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
            print(f"  SKIP {sid}: {size/1e9:.2f} GB — looks truncated (expect ~2.8 GB)")
            continue
        yield sid, path


def sample_weights(ts: np.ndarray) -> np.ndarray:
    dt = np.gradient(ts)
    return np.clip(dt, 0.0, MAX_GAP_S)


def in_spontaneous(ts: np.ndarray, starts: np.ndarray, stops: np.ndarray) -> np.ndarray:
    order = np.argsort(starts)
    starts, stops = starts[order], stops[order]
    assert np.all(starts[1:] >= stops[:-1] - 1e-9), "spontaneous intervals overlap"
    idx = np.searchsorted(starts, ts, side="right") - 1
    valid = idx >= 0
    mask = np.zeros(ts.shape, dtype=bool)
    mask[valid] = ts[valid] < stops[idx[valid]]
    return mask


def stationary_bouts(ts: np.ndarray, stat: np.ndarray, interval_id: np.ndarray) -> pd.DataFrame:

    flip = np.r_[True, (stat[1:] != stat[:-1]) | (interval_id[1:] != interval_id[:-1])]
    grp = np.cumsum(flip) - 1
    rows = []
    for g in np.unique(grp[stat]):
        sel = np.where(grp == g)[0]
        t0, t1 = ts[sel[0]], ts[sel[-1]]
        rows.append({
            "interval_id": int(interval_id[sel[0]]),
            "t_start": float(t0),
            "t_stop": float(t1),
            "duration_s": float(t1 - t0),
            "n_samples": int(len(sel)),
        })
    return pd.DataFrame(rows, columns=["interval_id", "t_start", "t_stop",
                                       "duration_s", "n_samples"])


def scan_session(sid: str, path: str) -> tuple[dict, pd.DataFrame]:
    with h5py.File(path, "r") as f:
        for p in (SPONT_START, SPONT_STOP, RUN_DATA, RUN_TIME):
            if p not in f:
                raise KeyError(f"{sid}: missing HDF5 path {p}")
        starts = np.asarray(f[SPONT_START][:], dtype=float)
        stops = np.asarray(f[SPONT_STOP][:], dtype=float)
        speed = np.asarray(f[RUN_DATA][:], dtype=float)
        ts = np.asarray(f[RUN_TIME][:], dtype=float)

    if len(speed) != len(ts):
        raise ValueError(f"{sid}: running data/timestamps length mismatch "
                         f"({len(speed)} vs {len(ts)})")

    order = np.argsort(ts)
    ts, speed = ts[order], speed[order]
    if USE_ABS_SPEED:
        speed = np.abs(speed)

    finite = np.isfinite(speed) & np.isfinite(ts)
    ts, speed = ts[finite], speed[finite]

    mask = in_spontaneous(ts, starts, stops)
    ts_s, sp_s = ts[mask], speed[mask]
    if len(ts_s) < 2:
        raise ValueError(f"{sid}: no running samples inside spontaneous intervals")


    o = np.argsort(starts)
    iid = np.searchsorted(starts[o], ts_s, side="right") - 1

    w = sample_weights(ts_s)
    total_s = float(w.sum())
    running = sp_s >= RUN_THR
    stat = sp_s <= STAT_THR
    mid = ~running & ~stat

    bouts = stationary_bouts(ts_s, stat, iid)
    bouts.insert(0, "session_id", sid)
    fittable = bouts[bouts["duration_s"] >= MIN_BOUT_S]

    row = {
        "session_id": sid,
        "n_spont_intervals": int(len(starts)),
        "spont_declared_s": float((stops - starts).sum()),
        "spont_covered_s": total_s,
        "n_samples": int(len(ts_s)),
        "median_speed_cms": float(np.median(sp_s)),
        "frac_running": float(w[running].sum() / total_s),
        "frac_stationary": float(w[stat].sum() / total_s),
        "frac_intermediate": float(w[mid].sum() / total_s),
        "stationary_s": float(w[stat].sum()),
        "n_stat_bouts": int(len(bouts)),
        "n_stat_bouts_fittable": int(len(fittable)),
        "stationary_fittable_s": float(fittable["duration_s"].sum()),
        "longest_stat_bout_s": float(bouts["duration_s"].max()) if len(bouts) else 0.0,
    }
    return row, bouts


def main() -> None:
    print(f"Locomotion scope diagnostic — cache: {CACHE_DIR}")
    print(f"RUN >= {RUN_THR} cm/s | STAT <= {STAT_THR} cm/s | "
          f"fittable bout >= {MIN_BOUT_S}s | abs(speed)={USE_ABS_SPEED}\n")

    rows, all_bouts = [], []
    for sid, path in session_files(CACHE_DIR):
        try:
            row, bouts = scan_session(sid, path)
        except (KeyError, ValueError, OSError) as e:
            print(f"  FAIL {sid}: {type(e).__name__}: {e}")
            continue
        rows.append(row)
        all_bouts.append(bouts)
        print(f"  {sid}: running {row['frac_running']:6.1%} | "
              f"stationary {row['frac_stationary']:6.1%} | "
              f"mid {row['frac_intermediate']:5.1%} | "
              f"stat bouts >= {MIN_BOUT_S}s: {row['n_stat_bouts_fittable']:4d} "
              f"({row['stationary_fittable_s']:7.1f}s)")

    if not rows:
        print("\nNo sessions scanned. Check CACHE_DIR and file sizes.")
        return

    summary = pd.DataFrame(rows)
    bouts_df = pd.concat(all_bouts, ignore_index=True) if all_bouts else pd.DataFrame()
    summary.to_csv("locomotion_scope_summary.csv", index=False)
    bouts_df.to_csv("locomotion_scope_bouts.csv", index=False)


    tw = summary["spont_covered_s"]
    pooled_run = float((summary["frac_running"] * tw).sum() / tw.sum())
    pooled_stat = float((summary["frac_stationary"] * tw).sum() / tw.sum())
    pooled_fittable = int(summary["n_stat_bouts_fittable"].sum())
    pooled_fittable_s = float(summary["stationary_fittable_s"].sum())

    print("\n" + "=" * 74)
    print(f"POOLED ({len(summary)} sessions, {tw.sum():.0f}s spontaneous analysed)")
    print(f"  time-weighted running fraction    : {pooled_run:.1%}")
    print(f"  time-weighted stationary fraction : {pooled_stat:.1%}")
    print(f"  stationary bouts >= {MIN_BOUT_S}s        : {pooled_fittable} "
          f"({pooled_fittable_s:.1f}s total)")
    print("=" * 74)

    if pooled_fittable >= N_BOUTS_MIN and pooled_fittable_s >= 60.0:
        print("VERDICT: pooled stationary data may support an m-scan re-run (step15).")
        print("  -> step15 is a NEW result: criterion 5 quality gate applies, and the")
        print("     bin-width sweep must be repeated. Do not report it as a free bonus.")
        print("  -> Per-region unit counts still gate it (>= 20 good units).")
    else:
        print("VERDICT: insufficient stationary data for a credible re-run.")
        print("  -> §2b's claim narrows to AWAKE RUNNING mouse cortex. Say so plainly")
        print("     in Methods 2.1 and Results; the subcritical ceiling still stands,")
        print("     but its scope is smaller than currently written.")
    print("\nCaveat that travels either way: the retained-vs-declared spontaneous")
    print("coverage ratio below 1 means running samples don't span every interval.")
    print(summary[["session_id", "spont_declared_s", "spont_covered_s"]].to_string(index=False))


if __name__ == "__main__":
    main()
