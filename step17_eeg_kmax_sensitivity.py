
from __future__ import annotations

import sys
import os
import glob
import time
import logging

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step3b_branching_ratio_v2 as s3

logging.getLogger("mrestimator").setLevel(logging.ERROR)

DATA_DIR: str = sys.argv[1] if len(sys.argv) > 1 else "."


KMAX_MS_SWEEP: tuple[int, ...] = (960, 1920, 3840, 7680)
A7_PRIMARY_MS: int = 3840
MIN_KMAX_BINS: int = 20

CENSOR_REASON = "no decay within window"


def rules_for(bin_ms: float) -> list[tuple[str, int]]:
    out = [("current_60bins", 60)]
    for kms in KMAX_MS_SWEEP:
        kb = int(round(kms / bin_ms))
        out.append((f"A7_{kms}ms", max(MIN_KMAX_BINS, kb)))
    return out


def main() -> None:
    datainfo = os.path.join(DATA_DIR, "datainfo.mat")
    if not os.path.exists(datainfo):
        sys.exit(f"datainfo.mat not found in {DATA_DIR!r}")
    cmap = s3.condition_map(datainfo)
    set_files = sorted(glob.glob(os.path.join(DATA_DIR, "**", "*.set"), recursive=True))
    if not set_files:
        sys.exit(f"no .set files under {DATA_DIR!r}")

    import mne
    print(f"EEG KMAX sensitivity (A7.6) — {len(set_files)} .set files")
    print("Estimator: step3b.estimate_m, imported. Only KMAX_BINS varies.\n")

    rows: list[dict] = []
    t0 = time.time()
    for fi, f in enumerate(set_files, 1):
        stem = s3.norm(os.path.basename(f)[:-4])
        meta = cmap.get(stem)
        if meta is None:
            continue
        subject = os.path.basename(f).split("-")[0]

        ep = mne.io.read_epochs_eeglab(f, verbose="ERROR")
        fs = ep.info["sfreq"]
        X = ep.get_data()
        acts = [s3.population_activity(X[e]) for e in range(len(X))]
        rate = float(np.mean([a.sum() for a in acts]))

        dt = s3.choose_dt(acts)
        bin_ms = dt * 1000.0 / fs
        epoch_bins = min(len(s3.rebin(a, dt)) for a in acts) if acts else 0

        line = [f"  [{fi}/{len(set_files)}] {stem[:26]:26} {meta['condition']:8} "
                f"dt={dt:2d} bin={bin_ms:5.1f}ms"]
        for label, kmax_bins in rules_for(bin_ms):
            if epoch_bins <= kmax_bins + 2:
                rows.append(dict(subject=subject, level=meta["level"],
                                 condition=meta["condition"], conc=meta["conc"],
                                 rule=label, kmax_bins=kmax_bins, dt_samples=dt,
                                 bin_ms=bin_ms, m_MR=np.nan, tau_ms=np.nan,
                                 mr_ok=False, reason="epoch shorter than KMAX+2 bins",
                                 events_per_epoch=rate, n_epochs=len(X)))
                continue
            s3.KMAX_BINS = kmax_bins
            r = s3.estimate_m(acts, fs)
            tau_bins = r["tau_bins"]
            rows.append(dict(
                subject=subject, level=meta["level"], condition=meta["condition"],
                conc=meta["conc"], rule=label, kmax_bins=kmax_bins,
                window_ms=kmax_bins * bin_ms, tau_cap_ms=0.8 * kmax_bins * bin_ms,
                dt_samples=dt, bin_ms=bin_ms, m_MR=r["m_MR"], tau_bins=tau_bins,
                tau_ms=r["tau_ms"],
                window_over_tau=(kmax_bins / tau_bins) if np.isfinite(tau_bins) and tau_bins > 0 else np.nan,
                mr_ok=r["mr_ok"], reason=r["reason"],
                events_per_epoch=rate, n_epochs=len(X)))
            if label in ("current_60bins", f"A7_{A7_PRIMARY_MS}ms"):
                line.append(f"{label.split('_')[0]}: m={r['m_MR']:.3f}")
        print(" | ".join(line) + f"  ({time.time()-t0:.0f}s)")

    s3.KMAX_BINS = 60

    df = pd.DataFrame(rows)
    df.to_csv("eeg_kmax_sensitivity.csv", index=False)


    summ = []
    for rule, g in df.groupby("rule"):
        ok = g[g["mr_ok"] == True]
        cens = int((g["reason"] == CENSOR_REASON).sum())
        summ.append(dict(
            rule=rule, n_cells=len(g), n_mr_ok=len(ok), mr_ok_frac=len(ok) / len(g),
            censored_frac=cens / len(g),
            median_window_ms=g["window_ms"].median() if "window_ms" in g else np.nan,
            max_m=ok["m_MR"].max() if len(ok) else np.nan,
            median_m=ok["m_MR"].median() if len(ok) else np.nan,
            median_tau_ms=ok["tau_ms"].median() if len(ok) else np.nan,
            max_tau_ms=ok["tau_ms"].max() if len(ok) else np.nan,
            median_win_over_tau=ok["window_over_tau"].median() if len(ok) else np.nan,
            min_win_over_tau=ok["window_over_tau"].min() if len(ok) else np.nan,
            n_ge_098=int((ok["m_MR"] >= 0.98).sum()) if len(ok) else 0))
    order = ["current_60bins"] + [f"A7_{k}ms" for k in KMAX_MS_SWEEP]
    s = pd.DataFrame(summ).set_index("rule").reindex([r for r in order if r in df["rule"].values]).reset_index()
    s.to_csv("eeg_kmax_summary.csv", index=False)

    print("\n" + "=" * 104)
    print("EEG ARM — DOES m MOVE WHEN THE WINDOW IS FIXED IN TIME?")
    print("=" * 104)
    print(s.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


    cur = s[s["rule"] == "current_60bins"]
    pri = s[s["rule"] == f"A7_{A7_PRIMARY_MS}ms"]
    print("\n" + "-" * 104)
    if len(cur) and len(pri) and np.isfinite(cur["median_m"].iloc[0]) and np.isfinite(pri["median_m"].iloc[0]):
        d = abs(float(pri["median_m"].iloc[0]) - float(cur["median_m"].iloc[0]))
        print(f"median m: current {float(cur['median_m'].iloc[0]):.4f} -> "
              f"A7_{A7_PRIMARY_MS}ms {float(pri['median_m'].iloc[0]):.4f}   |shift| = {d:.4f}")
        if d < 0.002:
            print("VERDICT: EEG arm CLEAN. A7.6 discharged. step12/13 stand as scored;")
            print("         §2a needs no change. Record the check; do not re-run step12/13.")
        else:
            print("VERDICT: EEG arm INHERITS THE BIAS. branching_summary_A7.csv is emitted")
            print("         below for a step12/13 RE-SCORE. The falsification verdict must be")
            print("         RE-STATED, not patched — criteria 1-3 were scored against biased m.")
            print("         The original numbers stay in the record with the defect named.")
    print("\nCONVERGENCE (judge on these, not on m):")
    print("  median_tau_ms should FLATTEN across A7 windows -> fit converged, trust it.")
    print("  median_tau_ms tracking KMAX_MS upward          -> chasing the window; neither")
    print("                                                    old nor new number is usable.")
    print("  min_win_over_tau >= 5 -> exponential adequately constrained. Below ~3 it is not.")
    print("  NOTE the current rule's window_ms is NOT constant across cells (dt is")
    print("  data-driven, 4-48 ms), so its median_window_ms hides a ~12x spread that")
    print("  correlates with event rate — and event rate correlates with condition.")
    print("-" * 104)


    a7 = df[df["rule"] == f"A7_{A7_PRIMARY_MS}ms"].copy()
    out = a7[["subject", "level", "condition", "conc", "m_MR", "tau_ms", "dt_samples",
              "mr_ok", "reason", "events_per_epoch", "n_epochs"]]
    out.to_csv("branching_summary_A7.csv", index=False)
    print(f"\nWrote branching_summary_A7.csv ({len(out)} cells, "
          f"{int(out['mr_ok'].sum())} mr_ok) — drop-in shape for step4/step12/step13.")
    print("Do NOT overwrite branching_summary_v2.csv. Both stay in the record.")

    bad = df[df["mr_ok"] == False]
    if len(bad):
        print("\nDrop reasons by rule (criterion 5 reporting):")
        print(bad.groupby(["rule", "reason"]).size().to_string())


if __name__ == "__main__":
    main()
