
import sys
import os
import glob
import time
import logging

import numpy as np
import pandas as pd
import h5py

logging.getLogger("mrestimator").setLevel(logging.ERROR)




CACHE_DIR: str = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/allen_cache")
BIN_MS: float = 4.0
KMAX_BINS: int = 60
MIN_TAU_BINS: float = 2.0
MAX_TAU_FRAC: float = 0.8
MIN_UNITS: int = 20
RUN_THRESH: float = 3.0
STAT_THRESH: float = 0.5
MIN_INTERVAL_S: float = 2.0

AMP_CUTOFF_MAX: float = 0.1
PRESENCE_MIN: float = 0.9
ISI_VIOL_MAX: float = 0.5






def load_unit_spike_trains(h: h5py.File) -> list[np.ndarray]:
    st = h["units/spike_times"][:]
    idx = h["units/spike_times_index"][:]
    starts = np.concatenate(([0], idx[:-1]))
    return [st[a:b] for a, b in zip(starts, idx)]


def unit_regions(h: h5py.File) -> np.ndarray:
    peak = h["units/peak_channel_id"][:]
    eid = h["general/extracellular_ephys/electrodes/id"][:]
    loc = h["general/extracellular_ephys/electrodes/location"][:]
    loc = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in loc])

    id_to_loc = dict(zip(eid.tolist(), loc.tolist()))
    return np.array([id_to_loc.get(int(p), "unknown") for p in peak])


def good_unit_mask(h: h5py.File) -> np.ndarray:
    q = h["units/quality"][:]
    q = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in q])
    amp = h["units/amplitude_cutoff"][:]
    pres = h["units/presence_ratio"][:]
    isi = h["units/isi_violations"][:]
    with np.errstate(invalid="ignore"):
        m = (q == "good") & (amp <= AMP_CUTOFF_MAX) & (pres >= PRESENCE_MIN) & (isi <= ISI_VIOL_MAX)
    return m


def spontaneous_windows(h: h5py.File) -> list[tuple[float, float]]:
    a = h["intervals/spontaneous_presentations/start_time"][:]
    b = h["intervals/spontaneous_presentations/stop_time"][:]
    return [(float(s), float(e)) for s, e in zip(a, b)]


def running_series(h: h5py.File) -> tuple[np.ndarray, np.ndarray]:
    t = h["processing/running/running_speed/timestamps"][:]
    v = np.abs(h["processing/running/running_speed/data"][:])
    return t, v


def split_windows_by_state(win: list[tuple[float, float]], t: np.ndarray,
                           v: np.ndarray) -> dict[str, list[tuple[float, float]]]:
    out: dict[str, list[tuple[float, float]]] = {"running": [], "stationary": []}
    for (a, b) in win:
        m = (t >= a) & (t <= b)
        if m.sum() < 2:
            continue
        tt, vv = t[m], v[m]
        for state in ("running", "stationary"):
            in_state = vv >= RUN_THRESH if state == "running" else vv <= STAT_THRESH
            idx = np.where(in_state)[0]
            if idx.size == 0:
                continue
            for s in np.split(idx, np.where(np.diff(idx) > 1)[0] + 1):
                lo, hi = tt[s[0]], tt[s[-1]]
                if hi - lo >= MIN_INTERVAL_S:
                    out[state].append((float(lo), float(hi)))
    return out


def population_activity(trains: list[np.ndarray], unit_idx: np.ndarray,
                        intervals: list[tuple[float, float]],
                        bin_s: float) -> list[np.ndarray]:
    acts: list[np.ndarray] = []
    for (a, b) in intervals:
        nb = int((b - a) / bin_s)
        if nb < KMAX_BINS + 3:
            continue
        edges = a + np.arange(nb + 1) * bin_s
        pop = np.zeros(nb, dtype=np.int64)
        for ui in unit_idx:
            st = trains[ui]
            seg = st[(st >= edges[0]) & (st < edges[-1])]
            if seg.size:
                pop += np.histogram(seg, bins=edges)[0]
        acts.append(pop)
    return acts






def estimate_m(activities: list[np.ndarray], bin_s: float) -> dict:
    import mrestimator as mre
    out = dict(m_MR=np.nan, tau_ms=np.nan, mr_ok=False, reason="", n_bins=0)
    usable = [a for a in activities if len(a) > KMAX_BINS + 2]
    if len(usable) < 3:
        out["reason"] = "too few usable intervals"
        return out
    L = min(len(a) for a in usable)
    arr = np.vstack([a[:L] for a in usable])
    out["n_bins"] = int(arr.shape[1])
    if arr.sum() < 20:
        out["reason"] = "too few events"
        return out
    try:
        rk = mre.coefficients(arr, steps=(1, KMAX_BINS), dt=1, method="trialseparated")
        ft = mre.fit(rk, fitfunc="exponential")
        m, tau = float(ft.mre), float(ft.tau)
    except Exception as e:
        out["reason"] = f"fit failed: {type(e).__name__}"
        return out
    out.update(m_MR=m, tau_ms=tau * bin_s * 1000)
    if not np.isfinite(m) or not np.isfinite(tau):
        out["reason"] = "non-finite"
    elif tau < MIN_TAU_BINS:
        out["reason"] = "degenerate (tau too small)"
    elif tau > MAX_TAU_FRAC * KMAX_BINS:
        out["reason"] = "no decay within window"
    elif not (0 < m < 1.5):
        out["reason"] = "m out of range"
    else:
        out["mr_ok"] = True
        out["reason"] = "ok"
    return out






def scan_session(path: str, bin_s: float) -> list[dict]:
    records: list[dict] = []
    with h5py.File(path, "r") as h:
        sid = os.path.basename(path).replace("session_", "").replace(".nwb", "")
        trains = load_unit_spike_trains(h)
        regions = unit_regions(h)
        good = good_unit_mask(h)






        ivs = spontaneous_windows(h)

        good_regions = {}
        for r in np.unique(regions[good]):
            if isinstance(r, str) and r not in ("unknown",):
                uidx = np.where((regions == r) & good)[0]
                if uidx.size >= MIN_UNITS:
                    good_regions[r] = uidx
        print(f"  session {sid}: {len(good_regions)} regions >= {MIN_UNITS} units, "
              f"{len(ivs)} spontaneous windows (pooled)")

        for region, uidx in good_regions.items():
            acts = population_activity(trains, uidx, ivs, bin_s)
            r = estimate_m(acts, bin_s)
            rate = float(np.mean([a.sum() for a in acts])) if acts else 0.0
            records.append(dict(session=sid, region=region,
                                n_units=int(uidx.size), m_MR=r["m_MR"],
                                tau_ms=r["tau_ms"], n_bins=r["n_bins"],
                                mr_ok=r["mr_ok"], reason=r["reason"],
                                events_per_interval=rate))
    return records


BIN_SWEEP_MS: tuple[float, ...] = (4.0, 8.0, 16.0)


def main() -> None:
    files = sorted(glob.glob(os.path.join(CACHE_DIR, "session_*.nwb")))
    if not files:
        sys.exit(f"no session_*.nwb files in {CACHE_DIR!r}")
    print(f"Bin-width sweep {BIN_SWEEP_MS} ms over {len(files)} session(s)\n")

    all_records: list[dict] = []
    t0 = time.time()
    for bin_ms in BIN_SWEEP_MS:
        bin_s = bin_ms / 1000.0
        print(f"--- bin = {bin_ms} ms ---")
        for f in files:
            recs = scan_session(f, bin_s)
            for r in recs:
                r["bin_ms"] = bin_ms
            all_records += recs
        print(f"    ({time.time()-t0:.0f}s elapsed)\n")

    df = pd.DataFrame.from_records(all_records)
    out = os.path.join(CACHE_DIR, "allen_m_scan_binsweep.csv")
    df.to_csv(out, index=False)
    print(f"Done in {time.time()-t0:.0f}s. {len(df)} cells total -> {out}\n")

    print("================ BIN-WIDTH SENSITIVITY ================")
    print("If the m-ceiling stays <=~0.97 across bins, the subcritical NO-GO is solid.")
    print("If a coarser bin lifts cortex to >=0.98, the ceiling was a bin artifact.\n")
    print(f"  {'bin(ms)':>7}  {'n_ok':>5}  {'median':>7}  {'max':>6}  {'n>=0.98':>8}  {'n>1.0':>6}")
    for bin_ms in BIN_SWEEP_MS:
        ok = df[(df["bin_ms"] == bin_ms) & (df["mr_ok"])]
        if len(ok):
            print(f"  {bin_ms:>7}  {len(ok):>5}  {ok.m_MR.median():>7.3f}  {ok.m_MR.max():>6.3f}"
                  f"  {int((ok.m_MR>=0.98).sum()):>8}  {int((ok.m_MR>1.0).sum()):>6}")
        else:
            print(f"  {bin_ms:>7}  {0:>5}  {'--':>7}  {'--':>6}  {'--':>8}  {'--':>6}")
    print("\nReminder: n=4 sessions, awake mouse spiking, pooled spontaneous. This is the")
    print("cross-modality replication arm -- read it against propofol EEG (max m ~0.96).")


if __name__ == "__main__":
    main()
