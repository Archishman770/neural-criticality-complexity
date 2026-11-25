from __future__ import annotations
import sys
import os
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings("ignore")

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
N_SHUFFLES = 10_000
ALPHA = 0.05
QUALITY_GATE_FRAC = 0.15
SEED = 20260715


FIXED_DT_RESULTS = "fixed_dt_test_results.csv"   # written by step20_fixed_dt_test.py


def score_fixed_dt(data_dir: str) -> tuple[pd.DataFrame, dict] | None:
    """Score criteria 1-3 under fixed binning (amendment A8), from step20's output.

    Under data-driven binning (step3b), dt ~ 1/event-rate enters both the complexity
    and criticality sides of the model. step20 re-estimates m at each fixed dt in
    5..11 samples and fits the same models, with 10,000 within-subject permutations
    per dt. This function only scores those results; it does not refit.

    Verdict rule (fixed before scoring, logged in docs/decision_log.md): a criterion
    PASSES only if it passes at every fixed dt, FAILS if it fails at every fixed dt,
    and is MIXED otherwise.

    Args:
        data_dir: folder holding fixed_dt_test_results.csv.

    Returns:
        (per-dt table, {criterion: verdict}) or None if step20 has not been run.
    """
    path = os.path.join(data_dir, FIXED_DT_RESULTS)
    if not os.path.exists(path):
        return None
    r = pd.read_csv(path)
    r = r[r.scheme.str.startswith("B.")].copy()
    r["dt"] = r.scheme.str.extract(r"dt = (\d+)", expand=False).astype(int)
    mp = r[r.measure == "mpr_c"].set_index("dt")
    lz = r[r.measure == "lz"].set_index("dt")
    t = pd.DataFrame({
        "mpr_beta": mp.beta, "mpr_p": mp.p, "mpr_p_perm": mp.p_perm,
        "mpr_beta_var": mp.beta_var, "mpr_p_var": mp.p_var,
        "mpr_atten_pct": mp.attenuation_pct,
        "lz_beta": lz.beta, "lz_p": lz.p, "lz_p_var": lz.p_var,
        "lz_atten_pct": lz.attenuation_pct,
        "mpr_within_r": mp.perm_r, "lz_within_r": lz.perm_r,
    })
    t["crit1_pass"] = (t.mpr_beta < 0) & (t.mpr_p < ALPHA) & (t.mpr_p_perm < ALPHA)
    t["crit2_pass"] = t.mpr_p_var < ALPHA
    # "Equally strong or stronger", judged on the scale-free within-subject |r|, with
    # the 20% comparability band used in step13. Raw |beta| is NOT used here: LZ's SD
    # is ~3.5x MPR_C's, so a raw comparison always favours LZ (manuscript s2.5).
    t["crit3_lz_stronger"] = t.lz_within_r >= 0.8 * t.mpr_within_r

    def verdict(col: pd.Series, yes: str, no: str) -> str:
        return yes if col.all() else (no if not col.any() else "MIXED")

    v = dict(crit1=verdict(t.crit1_pass, "PASS", "FAIL"),
             crit2=verdict(t.crit2_pass, "PASS", "FAIL"),
             crit3=verdict(t.crit3_lz_stronger, "LIMITATION", "as hypothesized"))
    return t, v


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    comp = pd.read_csv(os.path.join(DATA_DIR, "complexity_summary.csv"))
    var = pd.read_csv(os.path.join(DATA_DIR, "variance_summary.csv"))
    br_all = pd.read_csv(os.path.join(DATA_DIR, "branching_summary_v2.csv"))
    br = br_all[br_all["mr_ok"] == True].copy()

    compm = comp.groupby(["subject", "condition"])[["lz", "mpr_c"]].mean().reset_index()
    varm = var.groupby(["subject", "condition"])["variance"].mean().reset_index()
    df = (compm.merge(varm, on=["subject", "condition"])
               .merge(br[["subject", "condition", "m_MR", "events_per_epoch"]],
                      on=["subject", "condition"]))
    df["dist"] = (1.0 - df["m_MR"]).abs()
    return df, br_all


def fit_dist(df: pd.DataFrame, measure: str, covars: list[str] | None = None
             ) -> tuple[float, float]:
    rhs = "dist" + ("".join(f" + {c}" for c in covars) if covars else "")
    m = smf.mixedlm(f"{measure} ~ {rhs}", df, groups=df["subject"]).fit()
    return float(m.params["dist"]), float(m.pvalues["dist"])


def permute_dist(df: pd.DataFrame, measure: str, n: int = N_SHUFFLES,
                 seed: int = SEED) -> tuple[float, np.ndarray]:
    rng = np.random.default_rng(seed)
    obs, _ = fit_dist(df, measure)
    work = df.copy()
    nulls = np.empty(n)
    ok = 0
    for i in range(n):


        work["dist"] = (work.groupby("subject")["dist"]
                            .transform(lambda s: rng.permutation(s.values)))
        try:
            b, _ = fit_dist(work, measure)
        except Exception:
            b = np.nan
        nulls[i] = b
        ok += int(np.isfinite(b))
        if (i + 1) % 1000 == 0:
            print(f"      ... {i+1}/{n} shuffles", flush=True)
    nulls = nulls[np.isfinite(nulls)]

    p = (1.0 + np.sum(nulls <= obs)) / (1.0 + len(nulls))
    return float(p), nulls


def main() -> None:
    df, br_all = load()
    n_cells, n_subj = len(df), df.subject.nunique()
    n_total_cells = len(br_all)
    n_ok = int(br_all["mr_ok"].sum())
    ok_frac = n_ok / n_total_cells if n_total_cells else 0.0

    print(f"Quality-passing cells joined: {n_cells} across {n_subj} subjects")
    print(f"m_MR range: {df.m_MR.min():.3f} to {df.m_MR.max():.3f}\n")

    results = {}


    crit5_pass = ok_frac >= QUALITY_GATE_FRAC
    results["crit5"] = dict(ok_frac=ok_frac, n_ok=n_ok, n_total=n_total_cells,
                            passed=crit5_pass)
    print(f"=== CRITERION 5 (quality gate) ===")
    print(f"   mr_ok: {n_ok}/{n_total_cells} = {ok_frac:.1%}  "
          f"(gate = {QUALITY_GATE_FRAC:.0%}) -> {'PASS' if crit5_pass else 'FAIL (underpowered)'}\n")


    b1, p1 = fit_dist(df, "mpr_c")
    crit1_sign_ok = b1 < 0
    crit1_pass = crit1_sign_ok and p1 < ALPHA
    results["crit1"] = dict(beta=b1, p=p1, sign_ok=crit1_sign_ok, passed=crit1_pass)
    print("=== CRITERION 1 (primary: mpr_c ~ dist + (1|subject)) ===")
    print(f"   beta_dist = {b1:+.4f}   p = {p1:.4f}")
    print(f"   sign negative as hypothesized? {crit1_sign_ok}")
    print(f"   -> {'PASS' if crit1_pass else 'FAIL'}\n")


    b2, p2 = fit_dist(df, "mpr_c", covars=["variance"])
    crit2_pass = p2 < ALPHA
    results["crit2"] = dict(beta=b2, p=p2, passed=crit2_pass)
    print("=== CRITERION 2 (confound: mpr_c ~ dist + variance + (1|subject)) ===")
    print(f"   beta_dist = {b2:+.4f}   p_dist = {p2:.4f}")
    print(f"   -> {'PASS' if crit2_pass else 'FAIL (dist lost significance with variance in model)'}\n")


    b3, p3 = fit_dist(df, "lz")
    b3v, p3v = fit_dist(df, "lz", covars=["variance"])
    lz_stronger = abs(b3) >= abs(b1) or p3 <= p1
    results["crit3"] = dict(beta_lz=b3, p_lz=p3, beta_lz_var=b3v, p_lz_var=p3v,
                            lz_stronger=lz_stronger)
    print("=== CRITERION 3 (LZ contrast) ===")
    print(f"   lz ~ dist:            beta = {b3:+.4f}  p = {p3:.4f}")
    print(f"   lz ~ dist + variance: beta = {b3v:+.4f}  p = {p3v:.4f}")
    print(f"   LZ as strong or stronger than MPR_C? {lz_stronger}")
    print(f"   -> {'REPORT AS LIMITATION (LZ matches/beats MPR_C)' if lz_stronger else 'MPR_C leads, as hypothesized'}\n")


    print(f"=== PERMUTATION ({N_SHUFFLES:,} within-subject shuffles; hard rule 3) ===")
    print("   [this is the slow part; shuffle count is NOT reduced]")
    pp_mpr, null_mpr = permute_dist(df, "mpr_c")
    print(f"   mpr_c: beta_dist={b1:+.4f}  parametric p={p1:.4f}  permutation p={pp_mpr:.4f}")
    pp_lz, null_lz = permute_dist(df, "lz", seed=SEED + 1)
    print(f"   lz   : beta_dist={b3:+.4f}  parametric p={p3:.4f}  permutation p={pp_lz:.4f}\n")
    results["perm"] = dict(p_perm_mpr=pp_mpr, p_perm_lz=pp_lz)


    reasons = br_all.loc[br_all["mr_ok"] != True, "reason"].value_counts()


    fixed = score_fixed_dt(DATA_DIR)
    if fixed is None:
        raise SystemExit(f"{FIXED_DT_RESULTS} not found in {DATA_DIR}. Run "
                         "step19_fixed_dt_branching.py then step20_fixed_dt_test.py first: "
                         "the binning-controlled arm (amendment A8) must be scored alongside "
                         "the pre-registered arm, never omitted.")
    fx, fxv = fixed
    print("=== BINNING-CONTROLLED ARM (amendment A8: fixed dt = 5..11 samples) ===")
    print(fx[["mpr_beta", "mpr_p", "mpr_p_perm", "mpr_p_var", "mpr_atten_pct",
              "lz_p", "lz_p_var", "lz_atten_pct"]].round(4).to_string())
    print(f"   Criterion 1: {fxv['crit1']}  ({int(fx.crit1_pass.sum())}/{len(fx)} dt)")
    print(f"   Criterion 2: {fxv['crit2']}  ({int(fx.crit2_pass.sum())}/{len(fx)} dt)")
    print(f"   Criterion 3: {fxv['crit3']}  ({int(fx.crit3_lz_stronger.sum())}/{len(fx)} dt)\n")
    results["fixed_dt"] = fxv

    h1_supported = crit1_pass and crit2_pass and crit5_pass and pp_mpr < ALPHA
    h1_supported_fixed = fxv["crit1"] == "PASS" and fxv["crit2"] == "PASS" and crit5_pass
    print("=== VERDICT (propofol arm) ===")
    print("   Binning-controlled spec (A8) -- the defensible one, reported first:")
    print(f"     Criterion 1: {fxv['crit1']}   Criterion 2: {fxv['crit2']}   "
          f"Criterion 3: {fxv['crit3']}")
    print(f"     -> {'meets' if h1_supported_fixed else 'DOES NOT meet'} criteria 1-2,5")
    print("   Pre-registered spec (data-driven dt, confounded by event rate):")
    print(f"   Criterion 1 (primary):        {'PASS' if crit1_pass else 'FAIL'}")
    print(f"   Criterion 2 (variance):       {'PASS' if crit2_pass else 'FAIL'}")
    print(f"   Criterion 3 (LZ contrast):    {'LIMITATION' if lz_stronger else 'as hypothesized'}")
    print(f"   Criterion 4 (offset tol.):    NOT REACHED (amendment A5 -- peak not bracketed)")
    print(f"   Criterion 5 (quality gate):   {'PASS' if crit5_pass else 'FAIL'}")
    print(f"   -> Propofol arm {'meets' if h1_supported else 'DOES NOT meet'} criteria 1-3,5 "
          f"for a per-dataset result.")
    print("   NOTE: per research_design Section 10, a per-dataset pass is evidence TOWARD H1,")
    print("   never confirmation of it -- H1 requires >=2 independent-axis datasets with a")
    print("   pre-specified offset tolerance, which A5 records as not reachable in vivo.")
    print("   NOTE: pipeline was post-hoc tuned, so this remains non-confirmatory.")


    rep = [
        "# Falsification pass/fail — propofol arm (Chennu propofol sedation EEG)",
        "",
        "Generated by `step12_falsification_passfail.py`. Discharges the debt logged in",
        "`decision_log.md` ('Explicit written falsification pass/fail ... remains outstanding')",
        "and in `research_design.md` amendment A2.",
        "",
        f"Cells: {n_cells} quality-passing (of {n_total_cells}) across {n_subj} subjects. ",
        f"m_MR range {df.m_MR.min():.3f}–{df.m_MR.max():.3f}.",
        "",
        "## 1. Pass/fail per criterion",
        "",
        "### 1a. Binning-controlled specification (amendment A8) -- reported first",
        "",
        "The pre-registered pipeline bins at the data-driven mean inter-event interval, so",
        "dt ~ 1/event-rate sits on both sides of the model (LZ vs events_per_epoch r = +0.94).",
        "Holding dt fixed (step19/step20) removes that path. Verdict rule: PASS only if the",
        "criterion passes at every fixed dt, FAIL if it fails at every dt, else MIXED.",
        "",
        "| dt | MPR_C beta | p | p_perm | p (+variance) | attenuation | LZ p (+variance) |",
        "|---|---|---|---|---|---|---|",
    ]
    for dt, row in fx.iterrows():
        rep.append(f"| {dt} | {row.mpr_beta:+.4f} | {row.mpr_p:.4f} | {row.mpr_p_perm:.4f} | "
                   f"{row.mpr_p_var:.3f} | {row.mpr_atten_pct:.0f}% | {row.lz_p_var:.3f} |")
    rep += [
        "",
        f"- Criterion 1 (primary): **{fxv['crit1']}** "
        f"({int(fx.crit1_pass.sum())}/{len(fx)} dt)",
        f"- Criterion 2 (variance confound): **{fxv['crit2']}** "
        f"({int(fx.crit2_pass.sum())}/{len(fx)} dt)",
        f"- Criterion 3 (LZ contrast): **{fxv['crit3']}** "
        f"({int(fx.crit3_lz_stronger.sum())}/{len(fx)} dt)",
        "",
        "### 1b. Pre-registered specification (data-driven dt)",
        "",
        "| Criterion | Test | Result | Verdict |",
        "|---|---|---|---|",
        f"| 1. Primary | `mpr_c ~ dist + (1\\|subject)` | beta={b1:+.4f}, p={p1:.4f} | "
        f"**{'PASS' if crit1_pass else 'FAIL'}** |",
        f"| 2. Confound | `+ variance` | beta={b2:+.4f}, p={p2:.4f} | "
        f"**{'PASS' if crit2_pass else 'FAIL'}** |",
        f"| 3. LZ contrast | `lz ~ dist` | beta={b3:+.4f}, p={p3:.4f} | "
        f"{'**LIMITATION** — LZ matches/beats MPR_C' if lz_stronger else 'MPR_C leads'} |",
        "| 4. Offset tolerance | cross-dataset | NOT REACHED | see amendment A5 |",
        f"| 5. Quality gate | mr_ok fraction | {ok_frac:.1%} (gate {QUALITY_GATE_FRAC:.0%}) | "
        f"**{'PASS' if crit5_pass else 'FAIL'}** |",
        "",
        "## 2. Confound-control result",
        "",
        f"With raw variance in the model, `dist` beta={b2:+.4f}, p={p2:.4f} "
        f"({'survives' if crit2_pass else 'does NOT survive'}).",
        "Per CLAUDE.md Section 2, the LZ/MPR_C asymmetry is a reported result: LZ carries a",
        "cleaner criticality signal (variance contributes nothing), MPR_C carries a criticality",
        "signal AND an independent amplitude component. State this plainly.",
        "",
        "## 3. Permutation p alongside parametric",
        "",
        f"- `mpr_c`: parametric p={p1:.4f}, permutation p={pp_mpr:.4f} "
        f"({N_SHUFFLES:,} within-subject shuffles)",
        f"- `lz`: parametric p={p3:.4f}, permutation p={pp_lz:.4f}",
        "",
        "## 4. mr_ok drop rate and reasons",
        "",
        f"Passed: {n_ok}/{n_total_cells} ({ok_frac:.1%}). Dropped: {n_total_cells - n_ok}.",
        "",
        "| reason | n |",
        "|---|---|",
    ]
    for r, c in reasons.items():
        rep.append(f"| {r} | {c} |")
    rep += [
        "",
        "## Standing caveats (must travel with this result)",
        "",
        "- The pipeline was post-hoc tuned; this is **not confirmatory**. Confirmatory status",
        "  requires the pre-registered model re-run on held-out data.",
        "- A per-dataset pass is **evidence toward H1, not confirmation** (research_design",
        "  Section 10). H1 needs >=2 independent-axis datasets with a pre-specified offset",
        "  tolerance — recorded as not reachable in vivo (amendment A5).",
        "- Propofol EEG shows the **rising limb only** (m_MR tops out ~0.957); the complexity",
        "  peak is not bracketed, so the offset is not identifiable in this arm.",
        "",
    ]
    out_md = os.path.join(DATA_DIR, "step12_falsification_report.md")
    with open(out_md, "w") as f:
        f.write("\n".join(rep))

    pd.DataFrame([dict(
        criterion1_beta=b1, criterion1_p=p1, criterion1_pass=crit1_pass,
        criterion2_beta=b2, criterion2_p=p2, criterion2_pass=crit2_pass,
        criterion3_lz_beta=b3, criterion3_lz_p=p3, criterion3_lz_stronger=lz_stronger,
        criterion5_ok_frac=ok_frac, criterion5_pass=crit5_pass,
        perm_p_mpr=pp_mpr, perm_p_lz=pp_lz, n_cells=n_cells, n_subjects=n_subj,
        n_shuffles=N_SHUFFLES,
        fixed_dt_crit1=fxv["crit1"], fixed_dt_crit2=fxv["crit2"],
        fixed_dt_crit3=fxv["crit3"],
        fixed_dt_crit2_max_p=float(fx.mpr_p_var.max()),
        fixed_dt_crit2_min_p=float(fx.mpr_p_var.min()),
    )]).to_csv(os.path.join(DATA_DIR, "step12_falsification_results.csv"), index=False)

    print(f"\nReport -> {out_md}")
    print(f"Results -> {os.path.join(DATA_DIR, 'step12_falsification_results.csv')}")


if __name__ == "__main__":
    main()
