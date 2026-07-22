
import sys, os, glob, time
import numpy as np

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
LEVELS = {1: "baseline", 2: "mild", 3: "moderate", 4: "recovery"}


def norm(n): return n.replace(" ", "_").replace(".", "_")


def condition_map(p):
    from pymatreader import read_mat
    return {norm(r[0]): dict(level=int(r[1]), condition=LEVELS[int(r[1])], conc=r[2])
            for r in read_mat(p)["datainfo"]}


def main():
    import mne, pandas as pd
    di = os.path.join(DATA_DIR, "datainfo.mat")
    if not os.path.exists(di):
        sys.exit(f"datainfo.mat not found in {DATA_DIR!r}")
    cmap = condition_map(di)
    files = sorted(glob.glob(os.path.join(DATA_DIR, "**", "*.set"), recursive=True))
    print(f"Found {len(files)} files.\n")

    rows = []; t0 = time.time()
    for fi, f in enumerate(files, 1):
        stem = norm(os.path.basename(f)[:-4]); meta = cmap.get(stem)
        if meta is None:
            continue
        subject = os.path.basename(f).split("-")[0]
        ep = mne.io.read_epochs_eeglab(f, verbose="ERROR")
        X = ep.get_data()
        var = X.var(axis=2).mean(axis=0)
        for ci, ch in enumerate(ep.ch_names):
            rows.append(dict(subject=subject, level=meta["level"],
                             condition=meta["condition"], conc=meta["conc"],
                             channel=ch, variance=float(var[ci])))
        print(f"  [{fi}/{len(files)}] {stem[:30]:30} {meta['condition']:8} ({time.time()-t0:.0f}s)")

    df = pd.DataFrame.from_records(rows)
    out = os.path.join(DATA_DIR, "variance_summary.csv")
    df.to_csv(out, index=False)
    print(f"\nDone. {len(df)} rows -> {out}")


if __name__ == "__main__":
    main()
