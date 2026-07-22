
import sys, os
from typing import Sequence
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DATA_DIR: str = sys.argv[1] if len(sys.argv) > 1 else "."
N_PERM: int = 10_000
SEED: int = 20260702
MEASURES: tuple[str, ...] = ("lz", "mpr_c", "variance")


DIRECTION: dict[str, int] = {"lz": +1, "mpr_c": +1, "variance": 0}


def load(data_dir: str) -> pd.DataFrame:
    comp = pd.read_csv(os.path.join(data_dir, "complexity_summary.csv"))
    var = pd.read_csv(os.path.join(data_dir, "variance_summary.csv"))
    br = pd.read_csv(os.path.join(data_dir, "branching_summary_v2.csv"))
    br = br[br["mr_ok"] == True]

    compm = comp.groupby(["subject", "condition"])[["lz", "mpr_c"]].mean().reset_index()
    varm = var.groupby(["subject", "condition"])["variance"].mean().reset_index()
    df = (compm.merge(varm, on=["subject", "condition"])
                .merge(br[["subject", "condition", "m_MR"]], on=["subject", "condition"]))
    return df.reset_index(drop=True)


def mean_within_subject_rho(mmr: np.ndarray, vals: np.ndarray,
                            blocks: Sequence[np.ndarray]) -> tuple[float, int]:
    rhos: list[float] = []
    for b in blocks:
        x, y = mmr[b], vals[b]
        if x.size > 2 and np.unique(x).size > 1 and np.unique(y).size > 1:
            r = spearmanr(x, y).statistic
            if np.isfinite(r):
                rhos.append(float(r))
    return (float(np.mean(rhos)), len(rhos)) if rhos else (np.nan, 0)


def main() -> None:
    df = load(DATA_DIR)
    rng = np.random.default_rng(SEED)
    blocks = [np.asarray(idx, dtype=int) for _, idx in df.groupby("subject").groups.items()]
    mmr = df["m_MR"].to_numpy()

    print(f"Cells: {len(df)}   subjects: {len(blocks)}   perms: {N_PERM}")
    print("Statistic: mean within-subject Spearman(m_MR, measure)\n")

    rows = []
    for col in MEASURES:
        obs_vals = df[col].to_numpy()
        obs, nsub = mean_within_subject_rho(mmr, obs_vals, blocks)

        null = np.empty(N_PERM)
        for i in range(N_PERM):
            v = obs_vals.copy()
            for b in blocks:
                v[b] = rng.permutation(v[b])
            null[i], _ = mean_within_subject_rho(mmr, v, blocks)


        p_two = (np.sum(np.abs(null) >= abs(obs)) + 1) / (N_PERM + 1)
        p_pos = (np.sum(null >= obs) + 1) / (N_PERM + 1)
        d = DIRECTION[col]
        p_dir = p_pos if d > 0 else (np.nan if d == 0 else (np.sum(null <= obs) + 1) / (N_PERM + 1))

        rows.append(dict(measure=col, mean_rho=obs, n_subj=nsub,
                         null_mean=null.mean(), null_sd=null.std(),
                         p_two_sided=p_two,
                         p_directional=p_dir, direction=d))
        pdir_str = f"{p_dir:.4f}" if np.isfinite(p_dir) else "  n/a (control)"
        print(f"  {col:9}: rho={obs:+.3f} (n_subj={nsub})   "
              f"null={null.mean():+.3f}±{null.std():.3f}   "
              f"p_two={p_two:.4f}   p_dir={pdir_str}")

    out = os.path.join(DATA_DIR, "step6_permutation_results.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nNull distributions + p-values -> {out}")
    print("\nFALSIFICATION CHECK: NOT EVALUATED. research_design.md thresholds are")
    print("not loaded, so pass/fail against the pre-specified falsification criteria")
    print("cannot be stated. These p-values are the inputs to that check, not the check.")
    print("Reminder: pipeline was tuned post-hoc -> this is robustness, not confirmation.")


if __name__ == "__main__":
    main()
