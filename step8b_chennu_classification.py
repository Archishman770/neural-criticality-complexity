
from __future__ import annotations

import sys
import os

import numpy as np
import pandas as pd
from scipy.stats import beta




DATA_DIR: str = sys.argv[1] if len(sys.argv) > 1 else "."

N_TRIALS: int = 40
ALPHA: float = 0.05
Z: float = 1.959963984540054

LEVELS = {1: "baseline", 2: "mild", 3: "moderate", 4: "recovery"}
BASELINE_LVL: int = 1
MODERATE_LVL: int = 3


THRESHOLD: int = 20
THRESHOLD_ALT: int = 21

PUBLISHED_N_DROWSY: int = 7


def clopper_pearson(k: int, n: int, alpha: float = ALPHA) -> tuple[float, float]:
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def wilson(k: int, n: int, z: float = Z) -> tuple[float, float]:
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z / denom) * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return float(max(0.0, centre - half)), float(min(1.0, centre + half))


def wald(k: int, n: int, z: float = Z) -> tuple[float, float]:
    p = k / n
    half = z * np.sqrt(p * (1 - p) / n)
    return float(max(0.0, p - half)), float(min(1.0, p + half))


CI_METHODS = {"clopper_pearson": clopper_pearson, "wilson": wilson, "wald": wald}


def load_behaviour(data_dir: str) -> pd.DataFrame:
    from pymatreader import read_mat
    path = os.path.join(data_dir, "datainfo.mat")
    if not os.path.exists(path):
        sys.exit(f"datainfo.mat not found in {data_dir!r}")
    rows = read_mat(path)["datainfo"]
    recs = []
    for name, lvl, conc, rt, correct in rows:
        recs.append(dict(subject=str(name).split("-")[0], level=int(lvl),
                         correct=int(round(float(correct))), conc=float(conc)))
    df = pd.DataFrame(recs)
    bad = df[(df["correct"] < 0) | (df["correct"] > N_TRIALS)]
    if len(bad):
        sys.exit(f"correct-response counts outside [0, {N_TRIALS}]:\n{bad.to_string()}\n"
                 f"The paper specifies 40 trials. Reconcile before trusting anything here.")
    return df


def classify(base_k: int, mod_k: int, method: str) -> tuple[bool, tuple, tuple]:
    f = CI_METHODS[method]
    b_lo, b_hi = f(base_k, N_TRIALS)
    m_lo, m_hi = f(mod_k, N_TRIALS)
    drowsy = m_hi < b_lo
    return drowsy, (b_lo, b_hi), (m_lo, m_hi)


def main() -> None:
    df = load_behaviour(DATA_DIR)
    base = df[df["level"] == BASELINE_LVL].set_index("subject")["correct"]
    mod = df[df["level"] == MODERATE_LVL].set_index("subject")["correct"]
    subjects = sorted(set(base.index) & set(mod.index), key=lambda s: (len(s), s))
    if len(subjects) != 20:
        print(f"  WARNING: {len(subjects)} subjects with both baseline and moderate "
              f"(expected 20). Proceeding.")

    rows = []
    for s in subjects:
        bk, mk = int(base[s]), int(mod[s])
        rec = dict(subject=s, baseline_correct=bk, moderate_correct=mk,
                   baseline_hit=bk / N_TRIALS, moderate_hit=mk / N_TRIALS,
                   step8_unresponsive_le20=bool(mk <= THRESHOLD),
                   step8_unresponsive_le21=bool(mk <= THRESHOLD_ALT))
        for meth in CI_METHODS:
            drowsy, (b_lo, b_hi), (m_lo, m_hi) = classify(bk, mk, meth)
            rec[f"chennu_drowsy_{meth}"] = bool(drowsy)
            if meth == "clopper_pearson":
                rec.update(base_ci_lo=b_lo, base_ci_hi=b_hi,
                           mod_ci_lo=m_lo, mod_ci_hi=m_hi, ci_gap=b_lo - m_hi)
        rows.append(rec)

    out = pd.DataFrame(rows)
    out.to_csv("step8b_chennu_classification.csv", index=False)


    print("=" * 104)
    print("PER-SUBJECT CLASSIFICATION  (N_TRIALS = 40, 95% CI)")
    print("  Chennu rule: drowsy iff moderate CI is BELOW and non-overlapping with baseline CI")
    print("=" * 104)
    show = out[["subject", "baseline_correct", "moderate_correct", "base_ci_lo", "base_ci_hi",
                "mod_ci_lo", "mod_ci_hi", "ci_gap", "chennu_drowsy_clopper_pearson",
                "chennu_drowsy_wilson", "chennu_drowsy_wald",
                "step8_unresponsive_le20", "step8_unresponsive_le21"]]
    print(show.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


    print("\n" + "=" * 104)
    print("GROUP SIZES")
    print("=" * 104)
    print(f"  PUBLISHED (Chennu et al. 2016)      : {PUBLISHED_N_DROWSY} drowsy / "
          f"{len(subjects)-PUBLISHED_N_DROWSY} responsive")
    for meth in CI_METHODS:
        n_d = int(out[f"chennu_drowsy_{meth}"].sum())
        match = "  <-- matches published COUNT" if n_d == PUBLISHED_N_DROWSY else ""
        print(f"  reproduced, {meth:16s}: {n_d} drowsy / {len(subjects)-n_d} responsive{match}")
    n20 = int(out["step8_unresponsive_le20"].sum())
    n21 = int(out["step8_unresponsive_le21"].sum())
    print(f"  step8 primary  (correct <= {THRESHOLD})      : {n20} unresponsive / {len(subjects)-n20} responsive")
    print(f"  step8 variant  (correct <= {THRESHOLD_ALT})      : {n21} unresponsive / {len(subjects)-n21} responsive")
    print("\n  NOTE: matching COUNT is not matching SUBJECTS. See agreement below.")


    print("\n" + "=" * 104)
    print("SUBJECT-LEVEL AGREEMENT vs the pre-registered split")
    print("=" * 104)
    verdicts = {}
    for meth in CI_METHODS:
        ch = out[f"chennu_drowsy_{meth}"]
        for label, col in (("<=20 primary", "step8_unresponsive_le20"),
                           ("<=21 variant", "step8_unresponsive_le21")):
            st = out[col]
            disagree = out.loc[ch != st, "subject"].tolist()
            n_agree = int((ch == st).sum())
            kappa = _kappa(ch.to_numpy(), st.to_numpy())
            print(f"  {meth:16s} vs {label:13s}: {n_agree}/{len(out)} agree | "
                  f"kappa {kappa:5.3f} | differ on: {disagree if disagree else 'NONE'}")
            verdicts[(meth, label)] = len(disagree)


    cols = [f"chennu_drowsy_{m}" for m in CI_METHODS]
    stable = out[cols].nunique(axis=1).eq(1).all()
    print(f"\n  CI-method ambiguity: {'MOOT — all three methods agree on every subject.' if stable else 'LIVE — methods disagree; see per-subject table. Report the instability.'}")
    if not stable:
        wobble = out.loc[out[cols].nunique(axis=1) > 1, ["subject"] + cols]
        print(wobble.to_string(index=False))


    print("\n" + "-" * 104)
    key = ("clopper_pearson", "<=20 primary")
    if verdicts[key] == 0:
        print("VERDICT: the published classification IS the pre-registered primary split.")
        print("  -> Re-run NOTHING. Report the concordance in one Methods sentence:")
        print("     the absolute <=50% criterion recovers Chennu's within-subject binomial")
        print("     grouping exactly on this cohort. The A4 hook stands as scored.")
    elif verdicts[("clopper_pearson", "<=21 variant")] == 0:
        print("VERDICT: the published classification equals the <=21 SENSITIVITY VARIANT,")
        print("  not the <=20 primary. Report BOTH. The pre-registered primary stays primary;")
        print("  the variant gains external authority as it matches the source paper. State")
        print("  plainly that the two differ by one subject and report the hook under both.")
    else:
        print(f"VERDICT: the published classification DIFFERS from the primary on "
              f"{verdicts[key]} subject(s).")
        print("  -> step8 re-runs under the Chennu grouping as a THIRD reported grouping.")
        print("     It does NOT replace the pre-registered <=20 primary. Report all three.")
        print("     With n=6, expect the hook's strength to move. Whatever it does, report it —")
        print("     including if it weakens or dies. That is what pre-registration is FOR.")
    print("-" * 104)
    print("\nStop here. This script classifies; it does not re-run anything.")


def _kappa(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a.astype(bool), b.astype(bool)
    n = len(a)
    po = float((a == b).mean())
    pe = float(a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean()))
    return 1.0 if pe == 1.0 else (po - pe) / (1 - pe)


if __name__ == "__main__":
    main()
