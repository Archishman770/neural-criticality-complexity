
import sys
import os
import time
import logging
from typing import Optional

import numpy as np
import pandas as pd

logging.getLogger("mrestimator").setLevel(logging.ERROR)




CACHE_DIR: str = sys.argv[1] if len(sys.argv) > 1 else "./allen_cache"
N_SESSIONS: int = 1
BIN_MS: float = 4.0
KMAX_BINS: int = 60
MIN_TAU_BINS: float = 2.0
MAX_TAU_FRAC: float = 0.8
MIN_UNITS: int = 20
RUN_THRESH: float = 3.0
STAT_THRESH: float = 0.5
MIN_INTERVAL_S: float = 2.0

UNIT_FILTERS = dict(amplitude_cutoff_maximum=0.1,
                    presence_ratio_minimum=0.9,
                    isi_violations_maximum=0.5)






def load_allen_cache(cache_dir: str):
    from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
    manifest = os.path.join(cache_dir, "manifest.json")
    return EcephysProjectCache.from_warehouse(manifest=manifest)


def allen_session_ids(cache, n: int) -> list[int]:
    tbl = cache.get_session_table()
    fc = tbl[tbl["session_type"] == "functional_connectivity"]
    use = fc if len(fc) >= n else tbl
    return list(use.index[:n])


def allen_spontaneous_intervals(session) -> list[tuple[float, float]]:
    sp = session.stimulus_presentations
    spont = sp[sp["stimulus_name"] == "spontaneous"]
    return [(float(r.start_time), float(r.stop_time)) for r in spont.itertuples()]


def _running_speed_series(session) -> tuple[np.ndarray, np.ndarray]:
    rs = session.running_speed
    t = ((rs["start_time"].to_numpy() + rs["end_time"].to_numpy()) / 2.0)
    v = np.abs(rs["velocity"].to_numpy())
    return t, v


def split_intervals_by_state(session, intervals: list[tuple[float, float]]
                             ) -> dict[str, list[tuple[float, float]]]:
    t, v = _running_speed_series(session)
    out: dict[str, list[tuple[float, float]]] = {"running": [], "stationary": []}
    for (a, b) in intervals:
        m = (t >= a) & (t <= b)
        if m.sum() < 2:
            continue
        tt, vv = t[m], v[m]
        for state, lo, hi in (("running", RUN_THRESH, np.inf),
                              ("stationary", -np.inf, STAT_THRESH)):
            in_state = (vv >= lo) & (vv < hi) if state == "running" else (vv <= hi)

            idx = np.where(in_state)[0]
            if idx.size == 0:
                continue
            splits = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
            for s in splits:
                lo_t, hi_t = tt[s[0]], tt[s[-1]]
                if hi_t - lo_t >= MIN_INTERVAL_S:
                    out[state].append((float(lo_t), float(hi_t)))
    return out


def allen_population_activity(session, unit_ids: list[int],
                              intervals: list[tuple[float, float]],
                              bin_s: float) -> list[np.ndarray]:
    spikes = session.spike_times
    acts: list[np.ndarray] = []
    for (a, b) in intervals:
        nb = int((b - a) / bin_s)
        if nb < KMAX_BINS + 3:
            continue
        edges = a + np.arange(nb + 1) * bin_s
        pop = np.zeros(nb, dtype=np.int64)
        for uid in unit_ids:
            st = spikes.get(uid)
            if st is None or st.size == 0:
                continue
            seg = st[(st >= edges[0]) & (st < edges[-1])]
            if seg.size:
                pop += np.histogram(seg, bins=edges)[0]
        acts.append(pop)
    return acts


def estimate_m(activities: list[np.ndarray], bin_s: float) -> dict:
    import mrestimator as mre
    out = dict(m_MR=np.nan, tau_bins=np.nan, tau_ms=np.nan,
               mr_ok=False, reason="", n_bins=0)
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

    out.update(m_MR=m, tau_bins=tau, tau_ms=tau * bin_s * 1000)
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


def region_units(session) -> dict[str, list[int]]:
    u = session.units
    for col, thr in (("amplitude_cutoff", UNIT_FILTERS["amplitude_cutoff_maximum"]),
                     ("presence_ratio", UNIT_FILTERS["presence_ratio_minimum"]),
                     ("isi_violations", UNIT_FILTERS["isi_violations_maximum"])):
        if col in u.columns:
            u = u[u[col] <= thr] if col != "presence_ratio" else u[u[col] >= thr]
    out: dict[str, list[int]] = {}
    for region, grp in u.groupby("ecephys_structure_acronym"):
        if isinstance(region, str) and len(grp) >= MIN_UNITS:
            out[region] = list(grp.index)
    return out


def main() -> None:
    bin_s = BIN_MS / 1000.0
    cache = load_allen_cache(CACHE_DIR)
    sids = allen_session_ids(cache, N_SESSIONS)
    print(f"Scanning {len(sids)} sessions, bin={BIN_MS} ms\n")

    records = []
    t0 = time.time()
    for si, sid in enumerate(sids, 1):
        session = cache.get_session_data(sid)
        spont = allen_spontaneous_intervals(session)
        if not spont:
            print(f"  [{si}/{len(sids)}] session {sid}: no spontaneous epochs, skip")
            continue
        by_state = split_intervals_by_state(session, spont)
        runits = region_units(session)
        print(f"  [{si}/{len(sids)}] session {sid}: {len(runits)} regions, "
              f"run={len(by_state['running'])} stat={len(by_state['stationary'])} intervals")

        for region, uids in runits.items():
            for state, ivs in by_state.items():
                if len(ivs) < 3:
                    continue
                acts = allen_population_activity(session, uids, ivs, bin_s)
                r = estimate_m(acts, bin_s)
                rate = float(np.mean([a.sum() for a in acts])) if acts else 0.0
                records.append(dict(session=sid, region=region, state=state,
                                    n_units=len(uids), m_MR=r["m_MR"],
                                    tau_ms=r["tau_ms"], n_bins=r["n_bins"],
                                    mr_ok=r["mr_ok"], reason=r["reason"],
                                    events_per_interval=rate))
        print(f"      ({time.time()-t0:.0f}s elapsed)")

    df = pd.DataFrame.from_records(records)
    out = os.path.join(CACHE_DIR, "allen_m_scan.csv")
    df.to_csv(out, index=False)

    ok = df[df["mr_ok"]] if len(df) else df
    print(f"\nDone in {time.time()-t0:.0f}s. {len(df)} cells, {len(ok)} passed quality "
          f"-> {out}")
    if len(ok):
        print(f"\nm_MR (quality-passing): min={ok.m_MR.min():.3f}  "
              f"median={ok.m_MR.median():.3f}  max={ok.m_MR.max():.3f}")
        print("\n================ DECISION GATE ================")
        n_near = int((ok.m_MR >= 0.98).sum())
        n_super = int((ok.m_MR > 1.00).sum())
        print(f"  cells with m >= 0.98 (near-critical): {n_near} / {len(ok)}")
        print(f"  cells with m  > 1.00 (supercritical): {n_super} / {len(ok)}")
        if n_super >= 3 or n_near >= 5:
            print("  -> Criticality axis reaches the peak region. Offset is measurable")
            print("     here. GO: build the full spiking arm (add complexity + the")
            print("     falsification/variance/permutation stack from the EEG arm).")
        else:
            print("  -> Still topping out subcritical, like propofol. The awake mouse")
            print("     brain may not bracket the peak either. NO-GO on this dataset as")
            print("     the offset arm; escalate to a supercritical-capable preparation")
            print("     (seizure spiking, or Plenz-lineage slice/culture).")
        print("  (Reminder: this is the criticality axis only. Complexity, and the")
        print("   offset itself, are computed only if this gate says GO.)")
    else:
        print("\nNo quality-passing cells -- inspect the `reason` column before concluding.")


if __name__ == "__main__":
    main()
