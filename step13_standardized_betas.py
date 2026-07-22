from __future__ import annotations
import sys
import os
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings("ignore")

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "."


def load() -> pd.DataFrame:
    comp = pd.read_csv(os.path.join(DATA_DIR, "complexity_summary.csv"))
    var = pd.read_csv(os.path.join(DATA_DIR, "variance_summary.csv"))
    br = pd.read_csv(os.path.join(DATA_DIR, "branching_summary_v2.csv"))
    br = br[br["mr_ok"] == True].copy()

    compm = comp.groupby(["subject", "condition"])[["lz", "mpr_c"]].mean().reset_index()
    varm = var.groupby(["subject", "condition"])["variance"].mean().reset_index()
    df = (compm.merge(varm, on=["subject", "condition"])
               .merge(br[["subject", "condition", "m_MR"]], on=["subject", "condition"]))
    df["dist"] = (1.0 - df["m_MR"]).abs()
    return df


def main() -> None:
    df = load()
    print(f"Cells: {len(df)}  subjects: {df.subject.nunique()}\n")

    rows = []
    print(f"{'measure':>8} {'range':>8} {'raw_beta':>10} {'std_beta':>10} "
          f"{'p':>10} {'std_beta_+var':>14} {'p_+var':>10}")
    print("-" * 76)
    for meas in ["lz", "mpr_c"]:
        rng = float(df[meas].max() - df[meas].min())
        sd = float(df[meas].std())
        df[meas + "_z"] = (df[meas] - df[meas].mean()) / sd

        raw = smf.mixedlm(f"{meas} ~ dist", df, groups=df["subject"]).fit()
        std = smf.mixedlm(f"{meas}_z ~ dist", df, groups=df["subject"]).fit()
        stdv = smf.mixedlm(f"{meas}_z ~ dist + variance", df, groups=df["subject"]).fit()

        rows.append(dict(measure=meas, obs_range=rng, obs_sd=sd,
                         raw_beta=float(raw.params["dist"]),
                         raw_p=float(raw.pvalues["dist"]),
                         std_beta=float(std.params["dist"]),
                         std_p=float(std.pvalues["dist"]),
                         std_beta_with_variance=float(stdv.params["dist"]),
                         std_p_with_variance=float(stdv.pvalues["dist"])))
        print(f"{meas:>8} {rng:>8.3f} {raw.params['dist']:>10.4f} "
              f"{std.params['dist']:>10.4f} {std.pvalues['dist']:>10.2e} "
              f"{stdv.params['dist']:>14.4f} {stdv.pvalues['dist']:>10.2e}")

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA_DIR, "step13_standardized_betas.csv"), index=False)

    lz, mp = out.iloc[0], out.iloc[1]
    raw_ratio = abs(lz.raw_beta) / abs(mp.raw_beta)
    scale_ratio = lz.obs_sd / mp.obs_sd
    std_ratio = abs(lz.std_beta) / abs(mp.std_beta)

    print("\n--- READOUT ---")
    print(f"Raw beta ratio (LZ/MPR_C):        {raw_ratio:.2f}x")
    print(f"Measurement SD ratio (LZ/MPR_C):  {scale_ratio:.2f}x")
    print(f"STANDARDIZED beta ratio:          {std_ratio:.2f}x   <-- the meaningful one")
    print()
    if abs(std_ratio - 1.0) < 0.20:
        print("=> The raw-beta gap was a SCALE ARTIFACT. In SD units LZ and MPR_C track")
        print("   criticality COMPARABLY. Accurate sentence: 'LZ and MPR_C track")
        print("   criticality about equally well.' NOT 'LZ is 3x stronger.'")
        print("   Criterion 3 LIMITATION still stands (LZ is not weaker than MPR_C, so")
        print("   the structure-not-randomness claim is still undermined).")
    elif std_ratio > 1.20:
        print(f"=> LZ genuinely outperforms MPR_C by {std_ratio:.2f}x in SD units -- NOT a")
        print("   scale artifact. Criterion 3 LIMITATION stands, and stands strongly:")
        print("   the randomness measure beats the structure measure. Report plainly.")
    else:
        print(f"=> MPR_C outperforms LZ by {1/std_ratio:.2f}x in SD units. Re-check whether")
        print("   criterion 3's limitation flag should have triggered at all -- step12")
        print("   flagged it on raw betas, which this shows were not comparable.")

    print("\nNOTE: the variance-adjusted standardized betas above quantify the")
    print("amplitude-component asymmetry (CLAUDE.md Section 2) on a comparable scale:")
    print("compare how much each measure's beta shrinks when variance enters the model.")


if __name__ == "__main__":
    main()
