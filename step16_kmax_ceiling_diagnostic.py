
from __future__ import annotations

import sys
import os
import glob
import logging
import warnings
from typing import Iterator

import numpy as np
import pandas as pd
import h5py

logging.getLogger("mrestimator").setLevel(logging.CRITICAL)
warnings.filterwarnings("ignore")




CACHE_DIR: str = os.path.expanduser("~/allen_cache")


MAX_TAU_FRAC: float = 0.8
MIN_TAU_BINS: float = 2.0
M_CEIL_SCAN: float = 1.5
MIN_UNITS: int = 20
AMP_CUTOFF_MAX: float = 0.1
PRESENCE_MIN: float = 0.9
ISI_VIOL_MAX: float = 0.5
STAT_THR: float = 0.5
SMOOTH_WIN_S: float = 1.0


KMAX_SWEEP: tuple[int, ...] = (60, 120, 240, 480, 960)
T_TRIAL_S: float = 10.0


CENSOR_TOL: float = 0.05
CEIL_MARGIN: float = 0.005


SCANS: tuple[tuple[str, int], ...] = (
    ("pooled", 4), ("stationary", 4),
    ("pooled", 16), ("stationary", 16),
)

MAX_GAP_S: float = 0.1
SESSION_GLOB: str = "session_*.nwb"
EXPECTED_BYTES: int = 2_500_000_000

SPONT_START = "intervals/spontaneous_presentations/start_time"
SPONT_STOP = "intervals/spontaneous_presentations/stop_time"
RUN_DATA = "processing/running/running_speed/data"
RUN_TIME = "processing/running/running_speed/timestamps"
ELEC_ID = "general/extracellular_ephys/electrodes/id"
ELEC_LOC = "general/extracellular_ephys/electrodes/location"


def ceiling_for(kmax: int) -> float:
    return float(np.exp(-1.0 / (MAX_TAU_FRAC * kmax)))


def decode(a: np.ndarray) -> np.ndarray:
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def session_files(cache: str) -> Iterator[tuple[str, str]]:
    for p in sorted(glob.glob(os.path.join(cache, SESSION_GLOB))):
        sid = os.path.basename(p).replace("session_", "").replace(".nwb", "")
        if os.path.getsize(p) < EXPECTED_BYTES:
            print(f"  SKIP {sid}: truncated")
            continue
        yield sid, p


def load_units(f: h5py.File) -> pd.DataFrame:
    qual = decode(f["units/quality"][:])
    amp = np.asarray(f["units/amplitude_cutoff"][:], float)
    pres = np.asarray(f["units/presence_ratio"][:], float)
    isi = np.asarray(f["units/isi_violations"][:], float)
    pk = np.asarray(f["units/peak_channel_id"][:], np.int64)
    e_id = np.asarray(f[ELEC_ID][:], np.int64)
    e_loc = decode(f[ELEC_LOC][:])
    lmap = dict(zip(e_id.tolist(), e_loc.tolist()))
    region = np.array([lmap.get(int(c), "unknown") for c in pk])
    good = (qual == "good") & (amp <= AMP_CUTOFF_MAX) & \
           (pres >= PRESENCE_MIN) & (isi <= ISI_VIOL_MAX)
    return pd.DataFrame({"row": np.arange(len(qual)), "region": region, "good": good})


def unit_spikes(f: h5py.File) -> tuple[np.ndarray, np.ndarray]:
    st = np.asarray(f["units/spike_times"][:], float)
    idx = np.asarray(f["units/spike_times_index"][:], np.int64)
    return st, np.c_[np.r_[0, idx[:-1]], idx]


def in_spont(ts: np.ndarray, s0: np.ndarray, s1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    o = np.argsort(s0)
    a, b = s0[o], s1[o]
    idx = np.searchsorted(a, ts, side="right") - 1
    valid = idx >= 0
    mask = np.zeros(ts.shape, bool)
    mask[valid] = ts[valid] < b[idx[valid]]
    return mask, idx


def windows_for(path: str, partition: str) -> list[tuple[float, float]]:
    with h5py.File(path, "r") as f:
        s0 = np.asarray(f[SPONT_START][:], float)
        s1 = np.asarray(f[SPONT_STOP][:], float)
        if partition == "pooled":
            return list(zip(s0.tolist(), s1.tolist()))
        sp = np.asarray(f[RUN_DATA][:], float)
        ts = np.asarray(f[RUN_TIME][:], float)
    o = np.argsort(ts); ts, sp = ts[o], sp[o]
    ok = np.isfinite(ts) & np.isfinite(sp); ts, sp = ts[ok], sp[ok]
    fs = 1.0 / float(np.median(np.diff(ts)))
    n = max(1, int(round(SMOOTH_WIN_S * fs)))
    sm = np.convolve(sp, np.ones(n) / n, mode="same")
    mask, iid_all = in_spont(ts, s0, s1)
    ts_s, mag, iid = ts[mask], np.abs(sm)[mask], iid_all[mask]
    keep = mag <= STAT_THR
    if not keep.any():
        return []
    flip = np.r_[True, (keep[1:] != keep[:-1]) | (iid[1:] != iid[:-1])]
    grp = np.cumsum(flip) - 1
    out = []
    for g in np.unique(grp[keep]):
        sel = np.where(grp == g)[0]
        out.append((float(ts_s[sel[0]]), float(ts_s[sel[-1]])))
    return out


def segment(wins, t_trial: float) -> list[tuple[float, float]]:
    tr = []
    for t0, t1 in wins:
        for k in range(int((t1 - t0) // t_trial)):
            tr.append((t0 + k * t_trial, t0 + (k + 1) * t_trial))
    return tr


def fit_m(counts: np.ndarray, kmax: int) -> dict:
    import mrestimator as mre
    out = dict(m_MR=np.nan, tau_bins=np.nan, mr_ok=False, reason="")
    if counts.shape[0] < 3:
        out["reason"] = "too few trials"; return out
    if counts.sum() < 20:
        out["reason"] = "too few events"; return out
    try:
        rk = mre.coefficients(counts, steps=(1, kmax), dt=1, method="trialseparated")
        ft = mre.fit(rk, fitfunc="exponential")
        m, tau = float(ft.mre), float(ft.tau)
    except Exception as e:
        out["reason"] = f"fit failed: {type(e).__name__}"; return out
    out.update(m_MR=m, tau_bins=tau)
    if not np.isfinite(m) or not np.isfinite(tau):
        out["reason"] = "non-finite"
    elif tau < MIN_TAU_BINS:
        out["reason"] = "degenerate (tau too small)"
    elif tau > MAX_TAU_FRAC * kmax:
        out["reason"] = "no decay within window"
    elif not (0 < m < M_CEIL_SCAN):
        out["reason"] = "m out of physical range"
    else:
        out["mr_ok"] = True
    return out


def scan(sid: str, path: str, partition: str, bin_ms: int) -> list[dict]:
    trials = segment(windows_for(path, partition), T_TRIAL_S)
    if len(trials) < 3:
        return []
    bin_s = bin_ms / 1000.0
    n_bins = int(round(T_TRIAL_S / bin_s))
    rows = []
    with h5py.File(path, "r") as f:
        units = load_units(f)
        st, offs = unit_spikes(f)
        for region, g in units[units["good"]].groupby("region"):
            if region in ("unknown", "root", "") or len(g) < MIN_UNITS:
                continue
            spk = np.sort(np.concatenate([st[offs[r, 0]:offs[r, 1]]
                                          for r in g["row"].to_numpy()]))
            mat = np.empty((len(trials), n_bins), np.int64)
            for i, (t0, _) in enumerate(trials):
                edges = t0 + np.arange(n_bins + 1) * bin_s
                lo = np.searchsorted(spk, edges[0], "left")
                hi = np.searchsorted(spk, edges[-1], "left")
                mat[i] = np.histogram(spk[lo:hi], bins=edges)[0]
            for kmax in KMAX_SWEEP:
                if n_bins <= kmax + 2:
                    rows.append(dict(session_id=sid, region=region, partition=partition,
                                     bin_ms=bin_ms, kmax=kmax, m_MR=np.nan,
                                     tau_bins=np.nan, mr_ok=False,
                                     reason="trial shorter than KMAX+2 bins"))
                    continue
                r = fit_m(mat, kmax)
                rows.append(dict(
                    session_id=sid, region=region, partition=partition, bin_ms=bin_ms,
                    kmax=kmax, filter_ceiling=ceiling_for(kmax), n_units=int(len(g)),
                    n_trials=int(len(trials)),
                    tau_ms=r["tau_bins"] * bin_ms if np.isfinite(r["tau_bins"]) else np.nan,
                    mean_rate_hz=float(mat.sum() / (len(trials) * T_TRIAL_S)), **r))
    return rows


def main() -> None:
    cache = os.path.expanduser(sys.argv[1]) if len(sys.argv) > 1 else CACHE_DIR
    paths = list(session_files(cache))
    if not paths:
        print(f"No sessions in {cache}"); return

    print("KMAX ceiling diagnostic — the pre-registered rule (tau <= 0.8*KMAX) is HELD FIXED.")
    print("Induced m ceiling by KMAX:")
    for k in KMAX_SWEEP:
        print(f"    KMAX {k:4d} -> tau cap {MAX_TAU_FRAC*k:6.1f} bins -> m ceiling {ceiling_for(k):.5f}")
    print(f"\nPre-declared adequacy: censored < {CENSOR_TOL:.0%} AND max m <= ceiling - {CEIL_MARGIN}\n")

    rows = []
    for partition, bin_ms in SCANS:
        print(f"=== {partition} / {bin_ms} ms ===")
        for sid, path in paths:
            try:
                r = scan(sid, path, partition, bin_ms)
            except (KeyError, ValueError, OSError) as e:
                print(f"    FAIL {sid}: {type(e).__name__}: {e}"); continue
            rows.extend(r)
            print(f"    {sid}: {len(r)} cell-fits")

    if not rows:
        print("Nothing scanned."); return
    df = pd.DataFrame(rows)
    df.to_csv("allen_kmax_ceiling_diagnostic.csv", index=False)

    summ = []
    for (partition, bin_ms, kmax), g in df.groupby(["partition", "bin_ms", "kmax"]):
        ok = g[g["mr_ok"] == True]
        cens = int((g["reason"] == "no decay within window").sum())
        ceil = ceiling_for(kmax)
        adequate = (cens / len(g) < CENSOR_TOL) and \
                   (len(ok) > 0) and (ok["m_MR"].max() <= ceil - CEIL_MARGIN)
        summ.append(dict(
            partition=partition, bin_ms=bin_ms, kmax=kmax, filter_ceiling=ceil,
            n_cells=len(g), n_mr_ok=len(ok), mr_ok_frac=len(ok) / len(g),
            n_censored=cens, censored_frac=cens / len(g),
            max_m=ok["m_MR"].max() if len(ok) else np.nan,
            median_m=ok["m_MR"].median() if len(ok) else np.nan,
            gap_to_ceiling=(ceil - ok["m_MR"].max()) if len(ok) else np.nan,
            max_tau_ms=ok["tau_ms"].max() if len(ok) else np.nan,
            median_tau_ms=ok["tau_ms"].median() if len(ok) else np.nan,
            n_ge_098=int((ok["m_MR"] >= 0.98).sum()) if len(ok) else 0,
            n_gt_1=int((ok["m_MR"] > 1.0).sum()) if len(ok) else 0,
            adequate=adequate))
    s = pd.DataFrame(summ).sort_values(["partition", "bin_ms", "kmax"])
    s.to_csv("allen_kmax_ceiling_summary.csv", index=False)

    print("\n" + "=" * 100)
    print("KMAX SWEEP — does max m track the ceiling, or converge below it?")
    print("=" * 100)
    print(s[["partition", "bin_ms", "kmax", "filter_ceiling", "max_m", "gap_to_ceiling",
             "median_m", "censored_frac", "max_tau_ms", "median_tau_ms",
             "n_ge_098", "n_gt_1", "adequate"]].to_string(index=False,
             float_format=lambda x: f"{x:.4f}"))

    print("\n" + "-" * 100)
    print("READING THIS TABLE")
    print("  gap_to_ceiling ~ 0 at every KMAX      -> max m is PINNED. The ceiling is the")
    print("                                           filter. §2b's primary claim is an")
    print("                                           artifact; rewrite it.")
    print("  gap_to_ceiling GROWS, max_m flattens  -> real ceiling. §2b stands, stronger.")
    print("  max_tau_ms converges to a plausible cortical timescale (tens-hundreds of ms)")
    print("                                        -> fit is credible; trust the verdict.")
    print("  max_tau_ms grows in proportion to KMAX -> fit is chasing the window; NEITHER")
    print("                                           the old nor the new number is usable.")
    print("                                           Judge on tau_ms, not on m.")
    print("-" * 100)

    for (partition, bin_ms), g in s.groupby(["partition", "bin_ms"]):
        adq = g[g["adequate"]]
        if len(adq):
            k = int(adq["kmax"].min())
            print(f"  {partition}/{bin_ms}ms: smallest ADEQUATE KMAX = {k} "
                  f"(max m {float(adq[adq.kmax==k]['max_m'].iloc[0]):.4f}, "
                  f"ceiling {ceiling_for(k):.4f}) -> A7 candidate")
        else:
            print(f"  {partition}/{bin_ms}ms: NO adequate KMAX in sweep. Either the sweep")
            print(f"      is too narrow, or these timescales are unresolvable at "
                  f"T_TRIAL={T_TRIAL_S}s. That is itself the result — report it.")

    print("\nNEXT STEP IF THE FILTER WAS BINDING: write amendment A7 to research_design.md §9")
    print("BEFORE re-running the primary scan under the new KMAX. Log this diagnostic in")
    print("decision_log.md as exploratory. Do not fold the amendment into the re-run.")


if __name__ == "__main__":
    main()
