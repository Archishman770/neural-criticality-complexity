
import sys, os, glob, time, logging
import numpy as np
from scipy.signal import find_peaks

logging.getLogger("mrestimator").setLevel(logging.ERROR)




DATA_DIR    = sys.argv[1] if len(sys.argv) > 1 else "."
THRESHOLD   = 2.5
KMAX_BINS   = 60
MIN_TAU_BINS = 2.0
MAX_TAU_FRAC = 0.8
DT_BOUNDS   = (1, 12)
NEG_ONLY    = False
LEVELS      = {1: "baseline", 2: "mild", 3: "moderate", 4: "recovery"}


def norm(name):
    return name.replace(" ", "_").replace(".", "_")


def population_activity(epoch, thr=THRESHOLD, neg_only=NEG_ONLY):
    z = (epoch - epoch.mean(1, keepdims=True)) / (epoch.std(1, keepdims=True) + 1e-12)
    a = np.zeros(z.shape[1])
    for ch in range(z.shape[0]):
        sig = -z[ch] if neg_only else np.abs(z[ch])
        peaks, _ = find_peaks(sig, height=thr)
        a[peaks] += 1
    return a


def choose_dt(activities, bounds=DT_BOUNDS):
    ieis = []
    for a in activities:
        ev = np.where(a > 0)[0]
        if len(ev) > 1:
            ieis.append(np.diff(ev))
    if not ieis:
        return bounds[0]
    dt = int(round(np.mean(np.concatenate(ieis))))
    return int(np.clip(dt, *bounds))


def rebin(a, dt):
    nb = len(a) // dt
    return a[:nb * dt].reshape(nb, dt).sum(1) if nb else a[:0]


def estimate_m(activities, fs):
    import mrestimator as mre
    dt = choose_dt(activities)
    binned = [rebin(a, dt) for a in activities]
    binned = [b for b in binned if len(b) > KMAX_BINS + 2]
    out = dict(dt=dt, m_MR=np.nan, tau_bins=np.nan, tau_ms=np.nan,
               mr_ok=False, reason="")
    if len(binned) < 3:
        out["reason"] = "too few usable epochs"; return out
    L = min(len(b) for b in binned)
    arr = np.vstack([b[:L] for b in binned])
    if arr.sum() < 20:
        out["reason"] = "too few events"; return out
    try:
        rk = mre.coefficients(arr, steps=(1, KMAX_BINS), dt=1, method="trialseparated")
        ft = mre.fit(rk, fitfunc="exponential")
        m, tau = float(ft.mre), float(ft.tau)
    except Exception as e:
        out["reason"] = f"fit failed: {type(e).__name__}"; return out

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
        out["mr_ok"] = True; out["reason"] = "ok"
    return out


def condition_map(datainfo_path):
    from pymatreader import read_mat
    rows = read_mat(datainfo_path)["datainfo"]
    return {norm(r[0]): dict(level=int(r[1]), condition=LEVELS[int(r[1])],
                             conc=r[2]) for r in rows}


def main():
    import mne, pandas as pd
    datainfo = os.path.join(DATA_DIR, "datainfo.mat")
    if not os.path.exists(datainfo):
        sys.exit(f"datainfo.mat not found in {DATA_DIR!r}")
    cmap = condition_map(datainfo)
    set_files = sorted(glob.glob(os.path.join(DATA_DIR, "**", "*.set"), recursive=True))
    if not set_files:
        sys.exit(f"no .set files under {DATA_DIR!r}")
    print(f"Found {len(set_files)} .set files.\n")

    records = []; t0 = time.time()
    for fi, f in enumerate(set_files, 1):
        stem = norm(os.path.basename(f)[:-4])
        meta = cmap.get(stem)
        if meta is None:
            print(f"  [{fi}/{len(set_files)}] {stem}: not in datainfo, skipping"); continue
        subject = os.path.basename(f).split("-")[0]
        ep = mne.io.read_epochs_eeglab(f, verbose="ERROR")
        fs = ep.info["sfreq"]; X = ep.get_data()
        acts = [population_activity(X[e]) for e in range(len(X))]
        r = estimate_m(acts, fs)
        rate = float(np.mean([a.sum() for a in acts]))
        records.append(dict(subject=subject, level=meta["level"],
                            condition=meta["condition"], conc=meta["conc"],
                            m_MR=r["m_MR"], tau_ms=r["tau_ms"], dt_samples=r["dt"],
                            mr_ok=r["mr_ok"], reason=r["reason"],
                            events_per_epoch=rate, n_epochs=len(X)))
        flag = "ok " if r["mr_ok"] else "DROP"
        print(f"  [{fi}/{len(set_files)}] {stem[:30]:30} {meta['condition']:8} "
              f"m={r['m_MR']:.3f} tau={r['tau_ms']:.0f}ms [{flag}] ({time.time()-t0:.0f}s)")

    df = pd.DataFrame.from_records(records)
    out = os.path.join(DATA_DIR, "branching_summary_v2.csv")
    df.to_csv(out, index=False)
    n_ok = int(df["mr_ok"].sum())
    print(f"\nDone in {time.time()-t0:.0f}s.  {len(df)} cells, {n_ok} passed quality "
          f"({len(df)-n_ok} dropped) -> {out}")
    print("Next: step4 joins the mr_ok cells with complexity and tests the relationship,")
    print("controlling for event rate.")


if __name__ == "__main__":
    main()
