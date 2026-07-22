
import sys, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon
import statsmodels.formula.api as smf
import warnings; warnings.filterwarnings("ignore")

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "."


def zscore(s): return (s - s.mean()) / (s.std() + 1e-12)


def main():
    comp = pd.read_csv(os.path.join(DATA_DIR, "complexity_summary.csv"))
    var  = pd.read_csv(os.path.join(DATA_DIR, "variance_summary.csv"))
    br   = pd.read_csv(os.path.join(DATA_DIR, "branching_summary_v2.csv"))
    br   = br[br["mr_ok"] == True]

    compm = comp.groupby(["subject", "condition"])[["lz", "mpr_c"]].mean().reset_index()
    varm  = var.groupby(["subject", "condition"])["variance"].mean().reset_index()
    df = (compm.merge(varm, on=["subject", "condition"])
                .merge(br[["subject", "condition", "m_MR"]], on=["subject", "condition"]))
    df["dist"] = (1 - df["m_MR"]).abs()

    for c in ["lz", "mpr_c", "variance", "dist"]:
        df[c + "_z"] = df.groupby("subject")[c].transform(zscore)
    print(f"Cells: {len(df)}  subjects: {df.subject.nunique()}\n")

    print("=== A. complexity ~ dist + variance + (1|subject)  (does criticality survive?) ===")
    for meas in ["lz", "mpr_c"]:
        m = smf.mixedlm(f"{meas} ~ dist + variance", df, groups=df["subject"]).fit()
        print(f"   {meas:6}: dist p={m.pvalues['dist']:.3f}   variance p={m.pvalues['variance']:.3f}")

    print("\n=== B. variance ~ dist + (1|subject)  (does variance itself track criticality?) ===")
    m = smf.mixedlm("variance ~ dist", df, groups=df["subject"]).fit()
    print(f"   variance: beta_dist={m.params['dist']:+.4g}  p={m.pvalues['dist']:.3f}")

    print("\n=== C. within-subject: criticality vs complexity  vs  criticality vs variance ===")
    for label, col in [("complexity(mpr_c)", "mpr_c"), ("complexity(lz)", "lz"),
                       ("variance", "variance")]:
        rhos = [spearmanr(g["m_MR"], g[col])[0] for _, g in df.groupby("subject")
                if g["m_MR"].nunique() > 2]
        rhos = np.array([r for r in rhos if np.isfinite(r)])
        p = wilcoxon(rhos)[1] if len(rhos) >= 6 else np.nan
        print(f"   m_MR vs {label:18}: mean rho={rhos.mean():+.3f}  "
              f"{int((rhos>0).sum())}/{len(rhos)} pos  p={p:.3f}")

    print("\nRead: if dist dies in (A) but variance is significant, and variance tracks")
    print("criticality as well as complexity does in (C), the effect is an amplitude")
    print("confound -- not criticality-specific. If dist survives (A), it's real.")


if __name__ == "__main__":
    main()
