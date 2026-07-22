
import sys, os, glob, time
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from math import factorial


DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
D, TAU   = 6, 1
CHANNELS = None
LEVELS   = {1: "baseline", 2: "mild", 3: "moderate", 4: "recovery"}


def lz76(s):
    s = list(s); n = len(s); i = 0; c = l = k = k_max = 1
    while True:
        if s[i + k - 1] == s[l + k - 1]:
            k += 1
            if l + k > n:
                c += 1; break
        else:
            if k > k_max: k_max = k
            i += 1
            if i == l:
                c += 1; l += k_max
                if l + 1 > n: break
                i, k, k_max = 0, 1, 1
            else:
                k = 1
    return c

try:
    import antropy as _ant
    def lz_complexity(x):
        b = (np.asarray(x) > np.median(x)).astype(int)
        return _ant.lziv_complexity(b, normalize=True)
except ImportError:
    def lz_complexity(x):
        b = (np.asarray(x) > np.median(x)).astype(int)
        n = len(b)
        return lz76(b) * np.log2(n) / n

def _lehmer(P, d):
    facts = [factorial(d - 1 - i) for i in range(d)]
    codes = np.zeros(P.shape[0], dtype=np.int64)
    for i in range(d):
        codes += np.sum(P[:, i + 1:] < P[:, i:i + 1], axis=1) * facts[i]
    return codes

def mpr_complexity(x, d=D, tau=TAU):
    w = sliding_window_view(x, (d - 1) * tau + 1)[:, ::tau]
    P = np.argsort(w, axis=1, kind="stable")
    counts = np.bincount(_lehmer(P, d), minlength=factorial(d))
    N = len(counts); p = counts / counts.sum()
    nz = p > 0
    Sp = -np.sum(p[nz] * np.log(p[nz])); H = Sp / np.log(N)
    pe = 1.0 / N; m = (p + pe) / 2
    Sm = -np.sum(m * np.log(m)); Spe = np.log(N)
    Q = Sm - Sp / 2 - Spe / 2
    Q0 = -0.5 * ((N + 1) / N * np.log(N + 1) - 2 * np.log(2 * N) + np.log(N))
    return H, (Q / Q0) * H






def condition_map(datainfo_path):
    from pymatreader import read_mat
    rows = read_mat(datainfo_path)["datainfo"]
    out = {}
    for name, lvl, conc, rt, nc in rows:



        stem = name
        out[stem] = dict(level=int(lvl), condition=LEVELS[int(lvl)],
                         conc=conc, rt=rt, n_correct=nc)
    return out






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
        stem = os.path.basename(f)[:-4]
        meta = cmap.get(stem)
        if meta is None:
            print(f"  [{fi}/{len(set_files)}] {stem}: not in datainfo, skipping")
            continue
        subject = stem.split("-")[0]

        ep = mne.io.read_epochs_eeglab(f, verbose="ERROR")
        chans = CHANNELS or ep.ch_names
        data = ep.get_data(picks=chans)
        n_ep, n_ch, _ = data.shape

        for ci, ch in enumerate(chans):
            for ei in range(n_ep):
                sig = data[ei, ci]
                sig = sig - np.mean(sig)
                lz = lz_complexity(sig)
                H, C = mpr_complexity(sig)
                records.append(dict(subject=subject, file=stem,
                                    level=meta["level"], condition=meta["condition"],
                                    conc=meta["conc"], channel=ch, epoch=ei,
                                    lz=lz, H=H, mpr_c=C))
        dt = time.time() - t0
        print(f"  [{fi}/{len(set_files)}] {stem}  {meta['condition']:8} "
              f"{n_ep}ep x {n_ch}ch  ({dt:.0f}s elapsed)")

    df = pd.DataFrame.from_records(records)
    out_epochs = os.path.join(DATA_DIR, "complexity_epochs.csv")
    df.to_csv(out_epochs, index=False)

    summ = (df.groupby(["subject", "level", "condition", "conc", "channel"])
              [["lz", "H", "mpr_c"]].mean().reset_index())
    out_summary = os.path.join(DATA_DIR, "complexity_summary.csv")
    summ.to_csv(out_summary, index=False)

    print(f"\nDone in {time.time() - t0:.0f}s.")
    print(f"  {len(df):,} rows -> {out_epochs}")
    print(f"  {len(summ):,} rows -> {out_summary}")
    print("\nThis is the complexity axis. Next: branching ratio on the same files,")
    print("then join on (subject, condition) to test the offset hypothesis.")


if __name__ == "__main__":
    main()
