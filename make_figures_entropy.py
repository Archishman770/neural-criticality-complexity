"""Figures for the Entropy submission (methods framing). Every panel is built from a CSV
produced by a pipeline script; nothing is typed in by hand.

Fig 1  fit-window ceiling (Allen Neuropixels)      allen_kmax_ceiling_summary.csv (step16)
Fig 2  binning confound (propofol EEG)             Sedation-RestingState/*.csv (step3b, step20)
Fig 3  stationary MR fold (synthetic)              synthetic_mr_validation.csv (step9)
Fig 4  estimators under subsampling (synthetic)    transient_observed_validation.csv (step21),
                                                   avalanche_observed_validation.csv (step22),
                                                   kappa_rank_diagnostic_recordings.csv (step23)

Usage:  python3 make_figures_entropy.py      ->  figures_entropy/Fig{1..4}.{png,pdf}
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

import make_figures  # noqa: F401  (applies the shared rcParams)

EEG = "Sedation-RestingState"
OUT = "figures_entropy"
# Okabe-Ito subset, validated with the dataviz palette checker (all checks pass;
# orange is below 3:1 contrast, so every series also carries a distinct marker + legend).
BLUE, ORANGE, GREEN, VERM, GREY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#777777"


def save(fig: plt.Figure, name: str) -> None:
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{name}.png/.pdf")


def letter(ax: plt.Axes, s: str) -> None:
    ax.text(-0.2, 1.04, s, transform=ax.transAxes, fontsize=12, fontweight="bold")


def fig1() -> None:
    d = pd.read_csv("allen_kmax_ceiling_summary.csv")
    d = d[d.partition == "pooled"]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))
    ax = axes[0]
    k = np.geomspace(40, 1100, 200)
    ax.plot(k, np.exp(-1 / (0.8 * k)), color=GREY, lw=1.2, ls="--",
            label="filter ceiling exp(−1/(0.8·KMAX))")
    for b, col, mk in [(4, BLUE, "o"), (16, ORANGE, "s")]:
        s = d[(d.bin_ms == b) & d.max_m.notna()]
        ax.plot(s.kmax, s.max_m, color=col, marker=mk, ms=5, lw=1.5, label=f"max m, {b} ms bins")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("fit window KMAX (bins)")
    ax.set_ylabel("branching ratio m")
    ax.set_ylim(0.965, 1.0)
    ax.legend(fontsize=6.5, loc="lower right")
    letter(ax, "A")
    ax = axes[1]
    for b, col, mk in [(4, BLUE, "o"), (16, ORANGE, "s")]:
        s = d[(d.bin_ms == b) & (d.n_mr_ok > 0)]
        ax.plot(s.kmax, 100 * s.censored_frac, color=col, marker=mk, ms=5, lw=1.5, label=f"{b} ms bins")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("fit window KMAX (bins)")
    ax.set_ylabel("cells censored (%)")
    ax.legend(fontsize=7)
    letter(ax, "B")
    fig.tight_layout()
    save(fig, "Fig1_fit_window_ceiling")


def fig2() -> None:
    comp = (pd.read_csv(f"{EEG}/complexity_summary.csv")
            .groupby(["subject", "level"], as_index=False)[["lz", "mpr_c"]].mean())
    v2 = pd.read_csv(f"{EEG}/branching_summary_v2.csv")
    v2 = v2[v2.mr_ok == True].merge(comp, on=["subject", "level"])  # noqa: E712
    r = pd.read_csv(f"{EEG}/fixed_dt_test_results.csv")
    r["label"] = r.scheme.str.replace("A. data-driven dt (manuscript step3b)", "data-driven", regex=False)
    r["label"] = r.label.str.replace(r"B\. FIXED dt = (\d+) samples", r"fixed \1", regex=True)
    order = ["data-driven"] + [f"fixed {k}" for k in range(5, 12)]
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 2.8))
    ax = axes[0]
    rr = pearsonr(v2.events_per_epoch, v2.lz)[0]
    ax.scatter(v2.events_per_epoch, v2.lz, s=12, color=BLUE, alpha=0.75, edgecolors="none")
    ax.set_xlabel("events per epoch (event rate)")
    ax.set_ylabel("LZ complexity")
    ax.text(0.04, 0.92, f"r = {rr:+.2f}", transform=ax.transAxes, fontsize=8)
    letter(ax, "A")
    for ax, meas, name, let in [(axes[1], "mpr_c", "MPR$_C$", "B"), (axes[2], "lz", "LZ", "C")]:
        s = r[r.measure == meas].set_index("label").loc[order]
        x = np.arange(len(order))
        ax.plot(x, -np.log10(s.p), color=BLUE, marker="o", ms=5, lw=1.5, label="dist")
        ax.plot(x, -np.log10(s.p_var), color=VERM, marker="^", ms=5, lw=1.5, label="dist + variance")
        ax.axhline(-np.log10(0.05), color=GREY, lw=1, ls="--", label="p = 0.05")
        ax.set_xticks(x)
        ax.set_xticklabels(["data-\ndriven"] + [str(k) for k in range(5, 12)], fontsize=7)
        ax.set_xlabel("bin width Δt (samples)")
        ax.set_ylabel(f"−log$_{{10}}$ p, {name} ~ dist")
        ax.legend(fontsize=6.5, loc="upper right")
        letter(ax, let)
    fig.tight_layout()
    save(fig, "Fig2_binning_confound")


def fig3() -> None:
    d = pd.read_csv("synthetic_mr_validation.csv")
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    ax.plot([0.78, 1.32], [0.78, 1.32], color=GREY, lw=1, ls=":", label="identity")
    for f, col, mk, lab in [(1.0, BLUE, "o", "full observation"), (0.1, ORANGE, "s", "10% observation")]:
        s = d[d.subsample_f == f].sort_values("m_true")
        ax.plot(s.m_true, s.m_MR, color=col, marker=mk, ms=5, lw=1.5, label=lab)
    s = d[d.subsample_f == 1.0].sort_values("m_true")
    ax.plot(s.m_true, s.local_slope_pred, color=GREEN, lw=1.2, ls="--", label="local slope at fixed point")
    ax.axvline(1.0, color=GREY, lw=0.8)
    ax.set_xlabel("generative branching ratio m")
    ax.set_ylabel("stationary MR estimate")
    ax.legend(fontsize=6.5, loc="upper left")
    fig.tight_layout()
    save(fig, "Fig3_MR_fold")


def fig4() -> None:
    t = pd.read_csv("transient_observed_validation.csv")
    t = t[(t.world == "base") & (t.f == 0.05)]
    a = pd.read_csv("avalanche_observed_validation.csv")
    a = a[(a.world == "base") & (a.bin == 1)]
    k = pd.read_csv("kappa_rank_diagnostic_recordings.csv").dropna(subset=["kappa_4"])
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 2.8))
    ax = axes[0]
    for meth, thr, col, mk, lab in [("oracle", None, BLUE, "o", "onsets from full population"),
                                     ("observed", 2.0, VERM, "^", "onsets from observed units"),
                                     ("split", 2.0, ORANGE, "s", "split-half detection")]:
        s = t[(t.method == meth) & ((t.thr == thr) if thr else t.thr.isna())].sort_values("m_true")
        ax.plot(s.m_true, s.m_hat, color=col, marker=mk, ms=4.5, lw=1.5, label=lab)
    ax.axhline(1, color=GREY, lw=0.8)
    ax.axvline(1, color=GREY, lw=0.8)
    ax.set_xlabel("generative m")
    ax.set_ylabel("transient-onset estimate")
    ax.legend(fontsize=6, loc="lower right")
    letter(ax, "A")
    ax = axes[1]
    for f, col, mk in [(1.0, BLUE, "o"), (0.10, GREEN, "D"), (0.05, ORANGE, "s"), (0.02, VERM, "^")]:
        s = a[a.f == f].sort_values("m_true")
        ax.plot(s.m_true, s.kappa, color=col, marker=mk, ms=4.5, lw=1.5)
        ax.text(1.33, s.kappa.iloc[-1], f"f = {f:g}", fontsize=6.5, va="center")
    ax.axhline(1, color=GREY, lw=0.8)
    ax.axvline(1, color=GREY, lw=0.8)
    ax.set_xlabel("generative m")
    ax.set_ylabel("κ")
    ax.set_xlim(0.77, 1.42)
    ax.text(1.33, 1.17, "observed\nfraction", fontsize=6.5, va="center", color=GREY)
    letter(ax, "B")
    ax = axes[2]
    rho = spearmanr(k.kappa_4, k.mean_obs_rate)[0]
    ax.scatter(k.mean_obs_rate, k.kappa_4, s=8, color=BLUE, alpha=0.6, edgecolors="none")
    ax.set_xscale("log")
    ax.set_xlabel("mean observed count per 1 ms bin")
    ax.set_ylabel("κ (4 ms bins)")
    ax.text(0.04, 0.92, f"Spearman ρ = {rho:.3f}", transform=ax.transAxes, fontsize=7.5)
    letter(ax, "C")
    fig.tight_layout()
    save(fig, "Fig4_subsampled_estimators")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
