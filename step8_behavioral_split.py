
import sys
import os
from typing import Optional

import numpy as np
import pandas as pd




DATA_DIR: str = sys.argv[1] if len(sys.argv) > 1 else "."
# Primary rule is <= 20 (locked pre-inspection); pass 21 for the sensitivity variant.
THRESHOLD: int = int(sys.argv[2]) if len(sys.argv) > 2 else 20
THRESHOLD_ALT: int = 21
FIXED_DT_FILE: str = "branching_summary_fixeddt.csv"   # from step19 (amendment A8)
SEED_FIXED_DT: int = 20260927   # separate RNG stream so the original A4 numbers reproduce
N_PERM: int = 10_000
SEED: int = 20260704
MODERATE_LEVEL: int = 3

DIRECTION: dict[str, int] = {"mpr_c": +1, "lz": +1, "m_MR": +1, "variance": 0}






def load_behavior(datainfo_path: str, level: int) -> pd.DataFrame:
    from pymatreader import read_mat
    rows = read_mat(datainfo_path)["datainfo"]
    recs = []
    for name, lvl, conc, rt, correct in rows:
        if int(lvl) != level:
            continue
        subject = int(str(name).split("-")[0])
        recs.append(dict(subject=subject, correct=float(correct), conc=float(conc)))
    return pd.DataFrame.from_records(recs)


def build_table(data_dir: str, threshold: int) -> pd.DataFrame:
    comp = pd.read_csv(os.path.join(data_dir, "complexity_summary.csv"))
    var = pd.read_csv(os.path.join(data_dir, "variance_summary.csv"))
    br = pd.read_csv(os.path.join(data_dir, "branching_summary_v2.csv"))
    beh = load_behavior(os.path.join(data_dir, "datainfo.mat"), MODERATE_LEVEL)

    comp_m = comp[comp["level"] == MODERATE_LEVEL].groupby("subject")[["lz", "mpr_c"]].mean()
    var_m = var[var["level"] == MODERATE_LEVEL].groupby("subject")["variance"].mean()
    br_m = br[(br["level"] == MODERATE_LEVEL) & (br["mr_ok"] == True)].set_index("subject")["m_MR"]

    df = beh.set_index("subject").join(comp_m).join(var_m).join(br_m)
    df["responsive"] = (df["correct"] >= threshold + 1).astype(int)
    return df.reset_index()






def perm_diff(values: np.ndarray, resp: np.ndarray, direction: int,
              rng: np.random.Generator) -> dict:
    mask = np.isfinite(values)
    v, r = values[mask], resp[mask].astype(bool)
    n_r, n_u = int(r.sum()), int((~r).sum())
    if n_r < 2 or n_u < 2:
        return dict(n_resp=n_r, n_unresp=n_u, diff=np.nan, p_two=np.nan, p_dir=np.nan)
    obs = v[r].mean() - v[~r].mean()
    lab = r.copy()
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        rng.shuffle(lab)
        null[i] = v[lab].mean() - v[~lab].mean()
    p_two = (np.sum(np.abs(null) >= abs(obs)) + 1) / (N_PERM + 1)
    if direction > 0:
        p_dir = (np.sum(null >= obs) + 1) / (N_PERM + 1)
    elif direction < 0:
        p_dir = (np.sum(null <= obs) + 1) / (N_PERM + 1)
    else:
        p_dir = np.nan
    return dict(n_resp=n_r, n_unresp=n_u, diff=float(obs),
                p_two=float(p_two), p_dir=float(p_dir))


def cov_adjusted(df: pd.DataFrame, measure: str, rng: np.random.Generator) -> dict:
    import statsmodels.formula.api as smf
    sub = df[["responsive", "conc", measure]].dropna()
    if sub[measure].size < 6:
        return dict(reason="too few", coef=np.nan, p_perm=np.nan, conc_used=False)
    conc_varies = sub["conc"].nunique() > 1
    formula = f"{measure} ~ responsive + conc" if conc_varies else f"{measure} ~ responsive"
    obs = smf.ols(formula, sub).fit().params["responsive"]
    lab = sub["responsive"].to_numpy().copy()
    null = np.empty(N_PERM)
    work = sub.copy()
    for i in range(N_PERM):
        rng.shuffle(lab)
        work["responsive"] = lab
        null[i] = smf.ols(formula, work).fit().params["responsive"]
    p_perm = (np.sum(np.abs(null) >= abs(obs)) + 1) / (N_PERM + 1)
    return dict(reason="ok", coef=float(obs), p_perm=float(p_perm), conc_used=conc_varies)






def fixed_dt_m(data_dir: str, df: pd.DataFrame) -> pd.DataFrame | None:
    """Re-test the m_MR side of the split with fixed-dt branching ratios (amendment A8).

    The original m_MR uses step3b's data-driven bin width (dt ~ 1/event-rate), which
    A8 showed is confounded with the complexity measures. Here m is taken from
    step19 at each fixed dt = 5..11 samples and put through the same group permutation
    and concentration-adjusted permutation (10,000 shuffles each), on a separate RNG
    stream so the original results are unchanged.

    Args:
        data_dir: folder holding branching_summary_fixeddt.csv.
        df: per-subject moderate-condition table from build_table (responsive, conc).

    Returns:
        One row per dt, or None if step19 output is absent.
    """
    path = os.path.join(data_dir, FIXED_DT_FILE)
    if not os.path.exists(path):
        return None
    fx = pd.read_csv(path)
    fx = fx[(fx["level"] == MODERATE_LEVEL) & (fx["mr_ok"] == True)]  # noqa: E712
    rng = np.random.default_rng(SEED_FIXED_DT)
    rows = []
    for dt, g in fx.groupby("dt_fixed"):
        d = df.set_index("subject").join(g.set_index("subject")["m_MR"].rename("m_fixed"))
        d = d.reset_index()
        grp = perm_diff(d["m_fixed"].to_numpy(), d["responsive"].to_numpy(), +1, rng)
        adj = cov_adjusted(d, "m_fixed", rng)
        rows.append(dict(threshold=THRESHOLD, dt_fixed=int(dt), n_resp=grp["n_resp"],
                         n_unresp=grp["n_unresp"], diff=grp["diff"], p_two=grp["p_two"],
                         p_dir=grp["p_dir"], adj_coef=adj["coef"], adj_p_perm=adj["p_perm"]))
    return pd.DataFrame(rows)


def main() -> None:
    di = os.path.join(DATA_DIR, "datainfo.mat")
    if not os.path.exists(di):
        sys.exit(f"datainfo.mat not found in {DATA_DIR!r}")
    df = build_table(DATA_DIR, THRESHOLD)
    rng = np.random.default_rng(SEED)

    n_u = int((df["responsive"] == 0).sum())
    n_r = int((df["responsive"] == 1).sum())
    print("=" * 64)
    print("step8: behavioral split at MODERATE  (EXPLORATORY AMENDMENT)")
    print("=" * 64)
    print(f"Threshold: unresponsive = correct <= {THRESHOLD}")
    print(f"Split: {n_u} unresponsive / {n_r} responsive   (n=20 total)")
    print(f"Concentration at moderate: {sorted(df['conc'].dropna().unique())}")
    print(f"m_MR available (mr_ok) for {int(df['m_MR'].notna().sum())}/20 subjects\n")
    print("Directional hypothesis: unresponsive < responsive for mpr_c, lz, m_MR.")
    print("NOTE: n=6 in the unresponsive cell -> EXPLORATORY, not confirmatory.\n")

    print("--- Group comparison (permutation, 10k shuffles) ---")
    results = {}
    for meas in ["mpr_c", "lz", "m_MR", "variance"]:
        res = perm_diff(df[meas].to_numpy(), df["responsive"].to_numpy(),
                        DIRECTION[meas], rng)
        results[meas] = res
        d = res["diff"]
        if np.isfinite(d):
            pdir = f"{res['p_dir']:.4f}" if np.isfinite(res["p_dir"]) else "n/a"
            print(f"  {meas:9}: resp-unresp diff={d:+.4f}  "
                  f"p_two={res['p_two']:.4f}  p_dir={pdir}  "
                  f"(n {res['n_resp']}v{res['n_unresp']})")
        else:
            print(f"  {meas:9}: insufficient n ({res['n_resp']}v{res['n_unresp']})")

    print("\n--- Concentration-adjusted (does responsiveness survive controlling dose?) ---")
    cov = {}
    for meas in ["mpr_c", "lz", "m_MR"]:
        c = cov_adjusted(df, meas, rng)
        cov[meas] = c
        if c["reason"] == "ok":
            tag = "conc covariate" if c["conc_used"] else "conc constant, dropped"
            print(f"  {meas:9}: adj coef={c['coef']:+.4f}  p_perm={c['p_perm']:.4f}  [{tag}]")
        else:
            print(f"  {meas:9}: {c['reason']}")


    print("\n" + "=" * 64)
    print("SUCCESS CRITERION (pre-specified): MPR_C lower in unresponsive at")
    print("permutation p<.05 directional, surviving variance control + conc covariate.")
    print("=" * 64)
    mpr = results["mpr_c"]
    varr = results["variance"]
    covm = cov["mpr_c"]
    c1 = np.isfinite(mpr["p_dir"]) and mpr["p_dir"] < 0.05 and mpr["diff"] > 0
    c2 = covm["reason"] == "ok" and np.isfinite(covm["p_perm"]) and covm["p_perm"] < 0.05

    var_sep = np.isfinite(varr["p_two"]) and varr["p_two"] < 0.05
    print(f"  MPR_C directional p<.05 & correct sign : {'PASS' if c1 else 'FAIL'}")
    print(f"  Survives concentration covariate        : {'PASS' if c2 else 'FAIL'}")
    print(f"  Variance also separates groups?         : {'YES (weakens claim)' if var_sep else 'no'}")
    verdict = "SUPPORTED (exploratory)" if (c1 and c2 and not var_sep) else \
              "NOT SUPPORTED / underpowered"
    print(f"\n  AMENDMENT RESULT: {verdict}")
    print("  Reminder: n=6 unresponsive. Report as exploratory regardless of sign.")

    out = os.path.join(DATA_DIR, f"step8_behavioral_split_le{THRESHOLD}.csv")
    df.to_csv(out, index=False)
    print(f"\nPer-subject table -> {out}")

    # --- Amendment A8: m side of the dissociation under fixed binning -------
    fx = fixed_dt_m(DATA_DIR, df)
    if fx is None:
        print(f"\n[A8] {FIXED_DT_FILE} not found: the m_MR results above use the "
              "event-rate-confounded data-driven binning and must be reported with that caveat.")
        return
    print("\n" + "=" * 64)
    print("A8: m_MR responsive-vs-unresponsive under FIXED binning (dt = 5..11 samples)")
    print("=" * 64)
    print(fx[["dt_fixed", "n_resp", "n_unresp", "diff", "p_dir", "adj_coef",
              "adj_p_perm"]].round(4).to_string(index=False))
    n_grp = int((fx.p_dir < 0.05).sum())
    n_adj = int((fx.adj_p_perm < 0.05).sum())
    print(f"\n  group difference p_dir<.05 at {n_grp}/{len(fx)} dt; "
          f"survives concentration at {n_adj}/{len(fx)} dt")
    out_fx = os.path.join(DATA_DIR, f"step8_fixed_dt_le{THRESHOLD}.csv")
    fx.to_csv(out_fx, index=False)
    print(f"  -> {out_fx}")


if __name__ == "__main__":
    main()
