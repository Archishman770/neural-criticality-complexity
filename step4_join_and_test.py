
import sys, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon
import statsmodels.formula.api as smf
import warnings; warnings.filterwarnings("ignore")

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "."


def load():
    comp = pd.read_csv(os.path.join(DATA_DIR, "complexity_summary.csv"))
    br = pd.read_csv(os.path.join(DATA_DIR, "branching_summary_v2.csv"))
    br = br[br["mr_ok"] == True].copy()
    compm = comp.groupby(["subject", "condition"])[["lz", "H", "mpr_c"]].mean().reset_index()
    df = compm.merge(br[["subject", "condition", "m_MR", "events_per_epoch"]],
                     on=["subject", "condition"])
    df["dist"] = (1 - df["m_MR"]).abs()
    return df


def main():
    df = load()
    n_subj = df.subject.nunique()
    print(f"Quality-passing cells joined: {len(df)}  across {n_subj} subjects")
    print(f"m_MR range after cleaning: {df.m_MR.min():.3f} to {df.m_MR.max():.3f}\n")
    if len(df) < 12:
        print("Too few clean cells for a meaningful test. Inspect the drops in "
              "branching_summary_v2.csv (reason column).")
        return

    print("=== 1. MIXED-EFFECTS: complexity ~ dist + (1|subject) ===")
    for meas in ["lz", "mpr_c"]:
        m = smf.mixedlm(f"{meas} ~ dist", df, groups=df["subject"]).fit()
        print(f"   {meas:6}: beta_dist={m.params['dist']:+.4f}  p={m.pvalues['dist']:.3f}")

    print("\n=== 2. + event-rate control: complexity ~ dist + events_per_epoch + (1|subject) ===")
    for meas in ["lz", "mpr_c"]:
        m = smf.mixedlm(f"{meas} ~ dist + events_per_epoch", df, groups=df["subject"]).fit()
        print(f"   {meas:6}: beta_dist={m.params['dist']:+.4f}  p_dist={m.pvalues['dist']:.3f}"
              f"   (event-rate p={m.pvalues['events_per_epoch']:.3f})")

    print("\n=== 3. WITHIN-subject (does complexity follow m_MR across conditions?) ===")
    for meas in ["lz", "mpr_c"]:
        rhos = [spearmanr(g["m_MR"], g[meas])[0] for _, g in df.groupby("subject")
                if g["m_MR"].nunique() > 2]
        rhos = np.array([r for r in rhos if np.isfinite(r)])
        if len(rhos) >= 6:
            p = wilcoxon(rhos)[1]
            print(f"   {meas:6}: mean rho={rhos.mean():+.3f}  {int((rhos>0).sum())}/{len(rhos)} "
                  f"positive  Wilcoxon p={p:.3f}")

    print("\n=== 4. BETWEEN-subject (subject means) ===")
    sm = df.groupby("subject")[["lz", "mpr_c", "m_MR"]].mean()
    for meas in ["lz", "mpr_c"]:
        rho, p = spearmanr(sm["m_MR"], sm[meas])
        print(f"   {meas:6}: rho={rho:+.3f}  p={p:.3f}  (n={len(sm)})")


    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        cc = {"baseline": "#2c7fb8", "mild": "#7fcdbb", "moderate": "#fec44f", "recovery": "#de2d26"}
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
        for a, meas, lab in zip(ax, ["lz", "mpr_c"], ["Lempel-Ziv", "MPR statistical complexity"]):
            for c, col in cc.items():
                s = df[df.condition == c]
                a.scatter(s.m_MR, s[meas], c=col, s=42, alpha=.8, edgecolor="w", lw=.5, label=c)
            if len(df) > 2:
                z = np.polyfit(df.m_MR, df[meas], 1); xs = np.linspace(df.m_MR.min(), df.m_MR.max(), 50)
                a.plot(xs, np.polyval(z, xs), "k--", lw=1, alpha=.6)
            a.axvline(1, color="grey", ls=":", lw=1)
            a.set_xlabel("branching ratio m$_{MR}$ (1=critical)"); a.set_ylabel(meas); a.set_title(lab)
        ax[0].legend(fontsize=8, frameon=False)
        fig.suptitle(f"Quality-passing cells only (n={len(df)})")
        fig.tight_layout()
        out = os.path.join(DATA_DIR, "complexity_vs_criticality_clean.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"\nScatter -> {out}")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
