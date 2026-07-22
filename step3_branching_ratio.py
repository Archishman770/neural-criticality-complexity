
import sys, os, glob, time, logging
import numpy as np
from scipy.signal import find_peaks

logging.getLogger("mrestimator").setLevel(logging.ERROR)




DATA_DIR  = sys.argv[1] if len(sys.argv) > 1 else "."
THRESHOLD = 2.5
KMAX      = 50
NEG_ONLY  = False
LEVELS    = {1: "baseline", 2: "mild", 3: "moderate", 4: "recovery"}






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


def mr_branching_ratio(activities, kmax=KMAX):
    import mrestimator as mre
    arr = np.vstack(activities)
    if arr.sum() < 10:
        return np.nan, np.nan
    try:
        rk = mre.coefficients(arr, steps=(1, kmax), dt=1, method="trialseparated")
        ft = mre.fit(rk, fitfunc="exponential")
        return float(ft.mre), float(ft.tau)
    except Exception:
        return np.nan, np.nan


def naive_avalanche_sigma(activities):
    sigmas = []
    for a in activities:
        ev = np.where(a > 0)[0]
        if len(ev) < 4:
            continue
        dt = max(1, int(round(np.mean(np.diff(ev)))))
        nb = len(a) // dt
        if nb < 2:
            continue
        A = a[:nb * dt].reshape(nb, dt).sum(1)
        b = 0
        while b < len(A):
            if A[b] > 0:
                st = b
                while b < len(A) and A[b] > 0:
                    b += 1
                av = A[st:b]
                if len(av) >= 2 and av[0] > 0:
                    sigmas.append(av[1] / av[0])
            else:
                b += 1
    return float(np.mean(sigmas)) if sigmas else np.nan






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
        sys.exit(f"no .set files found under {DATA_DIR!r}")
    print(f"Found {len(set_files)} .set files.\n")

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
        X = ep.get_data()
        acts = [population_activity(X[e]) for e in range(len(X))]

        m, tau = mr_branching_ratio(acts)
        sigma = naive_avalanche_sigma(acts)
        rate = float(np.mean([a.sum() for a in acts]))

        records.append(dict(subject=subject, level=meta["level"],
                            condition=meta["condition"], conc=meta["conc"],
                            m_MR=m, tau_MR=tau, sigma_naive_BIASED=sigma,
                            events_per_epoch=rate, n_epochs=len(X)))
        print(f"  [{fi}/{len(set_files)}] {stem[:32]:32} {meta['condition']:8} "
              f"m_MR={m:.3f}  ({time.time()-t0:.0f}s)")

    df = pd.DataFrame.from_records(records)
    out = os.path.join(DATA_DIR, "branching_summary.csv")
    df.to_csv(out, index=False)
    print(f"\nDone in {time.time()-t0:.0f}s.  {len(df)} rows -> {out}")
    print("\nPrimary criticality measure = m_MR (branching ratio, subsampling-robust).")
    print("Join this with complexity_summary.csv on (subject, condition), then test")
    print("whether complexity peaks as m_MR approaches 1 -- your actual hypothesis.")


if __name__ == "__main__":
    main()
