"""Test the complexity-criticality relationship (manuscript s3.1) under FIXED binning.

Runs the manuscript's model specification against branching ratios estimated at a
constant bin width (from `step19_fixed_dt_branching.py`), rather than at the
data-driven bin width used in `step3b_branching_ratio_v2.py`.

Rationale: in step3b, dt = mean inter-event interval, so dt is essentially
1/event-rate. Both LZ and MPR_C also rise with event rate. The s3.1 relationship
therefore has the binning parameter on both sides. Holding dt fixed removes that
path; any surviving effect is attributable to the autocorrelation time tau, which
is the quantity of dynamical interest.

Applies the project's hard rules: the raw-variance confound control is reported
alongside every correlation, and cross-condition significance uses 10,000
within-subject label permutations.

Usage
-----
    python3 step20_fixed_dt_test.py <data_dir>
"""

from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
N_PERM = 10_000
RNG_SEED = 0
ORDER = ["baseline", "mild", "moderate", "recovery"]


def load_complexity(data_dir: str) -> pd.DataFrame:
    """Channel-averaged LZ / MPR_C / variance per subject-condition cell."""
    root = os.path.dirname(os.path.abspath(data_dir.rstrip("/")))
    comp_path = os.path.join(data_dir, "complexity_summary.csv")
    if not os.path.exists(comp_path):
        comp_path = os.path.join(root, "complexity_summary.csv")
    comp = pd.read_csv(comp_path)
    compm = comp.groupby(["subject", "level", "condition"], as_index=False)[
        ["lz", "H", "mpr_c"]
    ].mean()
    var = pd.read_csv(os.path.join(data_dir, "variance_summary.csv"))
    varm = var.groupby(["subject", "level"], as_index=False)["variance"].mean()
    return compm.merge(varm, on=["subject", "level"])


def within_subject_r(y: np.ndarray, x: np.ndarray, subj: np.ndarray) -> float:
    """Pearson r between subject-demeaned y and subject-demeaned x."""
    df = pd.DataFrame({"y": y, "x": x, "s": subj})
    yc = df.y - df.groupby("s").y.transform("mean")
    xc = df.x - df.groupby("s").x.transform("mean")
    return float(stats.pearsonr(yc, xc)[0])


def perm_test(d: pd.DataFrame, meas: str, n_perm: int = N_PERM) -> tuple[float, float]:
    """Within-subject label-permutation test on the complexity-dist coupling.

    Permutes condition labels within each subject (equivalently, shuffles the
    subject's dist values), preserving subject means and the complexity vector.

    Returns:
        (observed |r|, two-sided permutation p-value)
    """
    rng = np.random.default_rng(RNG_SEED)
    d = d.sort_values(["subject", "condition"])
    subj = d.subject.values
    obs = abs(within_subject_r(d[meas].values, d.dist.values, subj))

    # Pre-demean the complexity side once; permuting dist does not change it.
    ydf = pd.DataFrame({"y": d[meas].values, "s": subj})
    yc = (ydf.y - ydf.groupby("s").y.transform("mean")).values

    groups = [np.where(subj == s)[0] for s in pd.unique(subj)]
    dist = d.dist.values
    count = 0
    for _ in range(n_perm):
        p = dist.copy()
        for ix in groups:
            p[ix] = rng.permutation(p[ix])
        pdf = pd.DataFrame({"x": p, "s": subj})
        xc = (pdf.x - pdf.groupby("s").x.transform("mean")).values
        denom = np.sqrt((yc**2).sum() * (xc**2).sum())
        r = (yc * xc).sum() / denom if denom > 0 else 0.0
        if abs(r) >= obs:
            count += 1
    return obs, (count + 1) / (n_perm + 1)


def analyse(d: pd.DataFrame, label: str) -> list[dict]:
    """Fit the manuscript's models plus confound controls for one binning scheme."""
    rows = []
    print(f"\n{'='*78}\n{label}  (n={len(d)} cells, {d.subject.nunique()} subjects, "
          f"m range {d.m_MR.min():.3f}-{d.m_MR.max():.3f})\n{'='*78}")
    for meas in ["lz", "mpr_c"]:
        m1 = smf.mixedlm(f"{meas} ~ dist", d, groups=d["subject"]).fit()
        m2 = smf.mixedlm(f"{meas} ~ dist + variance", d, groups=d["subject"]).fit()
        b1, p1 = m1.params["dist"], m1.pvalues["dist"]
        b2, p2 = m2.params["dist"], m2.pvalues["dist"]
        obs_r, p_perm = perm_test(d, meas)
        atten = 100 * (1 - abs(b2 / b1)) if b1 != 0 else np.nan
        print(f"  {meas:6}  beta_dist = {b1:+.4f}  p = {p1:.3e}")
        print(f"          + variance  = {b2:+.4f}  p = {p2:.3e}   (attenuation {atten:+.0f}%)")
        print(f"          within-subj |r| = {obs_r:.3f}   p_perm({N_PERM}) = {p_perm:.5f}")
        rows.append(
            dict(scheme=label, measure=meas, beta=b1, p=p1, beta_var=b2, p_var=p2,
                 attenuation_pct=atten, perm_r=obs_r, p_perm=p_perm, n=len(d))
        )
    return rows


def main() -> None:
    comp = load_complexity(DATA_DIR)
    results = []

    # --- Reference: the manuscript's data-driven binning -------------------
    v2 = pd.read_csv(os.path.join(DATA_DIR, "branching_summary_v2.csv"))
    v2 = v2[v2.mr_ok == True].copy()  # noqa: E712
    d0 = comp.merge(v2[["subject", "level", "m_MR", "tau_ms", "dt_samples"]],
                    on=["subject", "level"])
    d0["dist"] = (1 - d0.m_MR).abs()
    results += analyse(d0, "A. data-driven dt (manuscript step3b)")

    # --- Fixed-dt sweep ----------------------------------------------------
    fx = pd.read_csv(os.path.join(DATA_DIR, "branching_summary_fixeddt.csv"))
    for dt in sorted(fx.dt_fixed.unique()):
        sub = fx[(fx.dt_fixed == dt) & (fx.mr_ok == True)]  # noqa: E712
        d = comp.merge(sub[["subject", "level", "m_MR", "tau_ms"]], on=["subject", "level"])
        d["dist"] = (1 - d.m_MR).abs()
        results += analyse(d, f"B. FIXED dt = {dt} samples")

    out = pd.DataFrame(results)
    path = os.path.join(DATA_DIR, "fixed_dt_test_results.csv")
    out.to_csv(path, index=False)

    print(f"\n\n{'='*78}\nSUMMARY: does the s3.1 effect survive fixed binning?\n{'='*78}")
    print(out.pivot_table(index="scheme", columns="measure",
                          values=["beta", "p", "p_perm"]).round(5).to_string())
    print(f"\n-> {path}")

    print(f"\n{'='*78}\nFALSIFICATION CRITERIA (fixed-dt arm)\n{'='*78}")
    fixed = out[out.scheme.str.startswith("B.")]
    for meas in ["lz", "mpr_c"]:
        s = fixed[fixed.measure == meas]
        n_sig = int((s.p < 0.05).sum())
        n_sig_var = int((s.p_var < 0.05).sum())
        n_neg = int(((s.p < 0.05) & (s.beta < 0)).sum())
        print(f"  {meas:6}: criterion 1 (beta<0, p<.05)     {n_neg}/{len(s)} fixed-dt values"
              f"  -> {'PASS' if n_neg > len(s)/2 else 'FAIL'}")
        print(f"          criterion 2 (survives variance) {n_sig_var}/{len(s)}"
              f"  -> {'PASS' if n_sig_var > len(s)/2 else 'FAIL'}")


if __name__ == "__main__":
    main()
