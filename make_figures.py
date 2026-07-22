from __future__ import annotations
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.titlesize": 10,
    "axes.labelsize": 9, "axes.linewidth": 0.8, "xtick.labelsize": 8,
    "ytick.labelsize": 8, "legend.fontsize": 8, "figure.dpi": 300,
    "savefig.dpi": 300, "savefig.bbox": "tight", "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False,
})
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00",
         grey="#999999", black="#000000", purple="#CC79A7", sky="#56B4E9")
PROJ = "."
OUT = "."


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}")
    plt.close(fig)
    print(f"wrote {name}.png / .pdf")


def fig1():
    bran = pd.read_csv(f"{PROJ}/branching_summary_v2.csv")
    comp = pd.read_csv(f"{PROJ}/complexity_summary.csv")
    bran = bran[bran.mr_ok == True].copy()
    cmean = comp.groupby(["subject", "condition"])[["lz", "mpr_c"]].mean().reset_index()
    df = bran.merge(cmean, on=["subject", "condition"], how="inner")
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.9))
    for ax, col, lab, colr in ((axes[0], "mpr_c", "MPR$_C$ (statistical complexity)", C["blue"]),
                               (axes[1], "lz", "LZ (algorithmic complexity)", C["orange"])):
        ax.scatter(df.m_MR, df[col], s=14, c=colr, alpha=0.55, edgecolors="none")
        z = np.polyfit(df.m_MR, df[col], 1)
        xs = np.linspace(df.m_MR.min(), df.m_MR.max(), 50)
        ax.plot(xs, np.polyval(z, xs), color=C["black"], lw=1.2, ls="--")
        ax.axvline(1.0, color=C["grey"], lw=0.8, ls=":")
        ax.set_xlabel("Branching ratio $m$")
        ax.set_ylabel(lab)
        ax.set_xlim(0.75, 1.02)
    axes[0].text(-0.18, 1.02, "A", transform=axes[0].transAxes, fontsize=13, fontweight="bold")
    axes[1].text(-0.18, 1.02, "B", transform=axes[1].transAxes, fontsize=13, fontweight="bold")
    fig.suptitle("Complexity rises as $m$ approaches criticality (rising limb only), propofol EEG",
                 fontsize=9.5, y=1.02)
    fig.tight_layout()
    save(fig, "Figure1_rising_limb")


def fig2():
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.0))
    ax = axes[0]
    tau = np.linspace(2, 200, 400)
    ax.plot(tau, np.exp(-1 / tau), color=C["blue"], lw=1.6)
    for kmax, colr, lab in [(48, C["red"], "KMAX=60 bins"), (768, C["green"], "KMAX=960 bins")]:
        ceil = np.exp(-1 / kmax)
        ax.axhline(ceil, color=colr, lw=1.0, ls="--")
        ax.axvline(kmax, color=colr, lw=0.8, ls=":")
        ax.text(200, ceil, f"  ceiling {ceil:.3f}\n  ({lab})", fontsize=6.5, va="center", color=colr)
    ax.set_xlabel(r"Fitted autocorrelation time $\tau$ (bins)")
    ax.set_ylabel(r"Estimated $m = \exp(-1/\tau)$")
    ax.set_title("The identity coupling filter to $m$", fontsize=8.5)
    ax.set_xlim(0, 200); ax.set_ylim(0.5, 1.0)
    ax.text(-0.2, 1.02, "A", transform=ax.transAxes, fontsize=13, fontweight="bold")

    ax = axes[1]
    groups = ["60 bins\n(original)", "3840 ms\n(A7), 4 ms", "3840 ms\n(A7), 16 ms"]









    def _panelB_from_csv():
        import os
        diag_p = f"{PROJ}/allen_kmax_ceiling_diagnostic.csv"
        scan_p = f"{PROJ}/allen_m_scan_stationary.csv"
        if not (os.path.exists(diag_p) and os.path.exists(scan_p)):
            return None
        diag = pd.read_csv(diag_p)
        scan = pd.read_csv(scan_p)

        d0 = diag[(diag.kmax == 60) & (diag.bin_ms == 4) & (diag.partition == "pooled")]
        ok0 = d0[d0.mr_ok == True]
        cens0 = 100.0 * (d0.reason == "no decay within window").mean()

        s = scan[(scan.run == "primary") & (scan.mr_ok == True)]
        s4, s16 = s[s.bin_ms == 4], s[s.bin_ms == 16]
        return dict(
            max_m=[float(ok0.m_MR.max()), float(s4.m_MR.max()), float(s16.m_MR.max())],
            med_m=[float(ok0.m_MR.median()), float(s4.m_MR.median()), float(s16.m_MR.median())],
            cens=[float(cens0), 0.0, 0.0],
            source="computed from allen_*.csv")

    vals = _panelB_from_csv()
    if vals is None:


        vals = dict(max_m=[0.979, 0.993, 0.975], med_m=[0.945, 0.978, 0.939],
                    cens=[28.3, 0.0, 0.0], source="FALLBACK: transcribed from Table 2")
    print(f"  Fig2 panel B source: {vals['source']}")
    max_m, med_m, cens = vals["max_m"], vals["med_m"], vals["cens"]

    x = np.arange(3); w = 0.35
    ax.bar(x - w / 2, max_m, w, color=C["blue"], label="max $m$")
    ax.bar(x + w / 2, med_m, w, color=C["sky"], label="median $m$")
    ax.axhline(np.exp(-1 / 48), color=C["red"], lw=1.0, ls="--")
    ax.text(2.4, np.exp(-1 / 48) + 0.001, "filter ceiling 0.979", color=C["red"],
            fontsize=6.5, ha="right", va="bottom")
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=7)
    ax.set_ylabel("Branching ratio $m$"); ax.set_ylim(0.90, 1.0)
    ax.legend(fontsize=7, loc="lower right")
    for xi, c in zip(x, cens):
        ax.text(xi, 0.902, f"{c:.0f}% censored", fontsize=6, ha="center",
                color=C["red"] if c > 0 else C["grey"])
    ax.set_title("Bin-specified filter caps $m$; time-specified does not", fontsize=8.5)
    ax.text(-0.2, 1.02, "B", transform=ax.transAxes, fontsize=13, fontweight="bold")
    fig.suptitle("The apparent subcritical ceiling was imposed by the fit-quality filter",
                 fontsize=9.5, y=1.03)
    fig.tight_layout()
    save(fig, "Figure2_fit_window_artifact")


def fig3():
    d = pd.read_csv(f"{PROJ}/synthetic_mr_validation.csv")
    full = d[d.subsample_f == 1.0].sort_values("m_true")
    sub = d[d.subsample_f == 0.1].sort_values("m_true")
    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    ax.plot([0.75, 1.35], [0.75, 1.35], color=C["grey"], lw=0.8, ls=":", label="identity")
    ax.axvline(1.0, color=C["grey"], lw=0.7, ls=":")
    ax.plot(full.m_true, full.m_MR, "o-", color=C["blue"], lw=1.4, ms=5, label="MR, full coverage")
    ax.plot(sub.m_true, sub.m_MR, "s--", color=C["sky"], lw=1.2, ms=4, label="MR, 10% coverage")
    ax.annotate("folds back\n(non-injective)", xy=(1.3, full[full.m_true == 1.3].m_MR.iloc[0]),
                xytext=(1.12, 0.72), fontsize=7,
                arrowprops=dict(arrowstyle="->", color=C["red"], lw=0.9), color=C["red"])
    ax.set_xlabel("True generative $m$"); ax.set_ylabel(r"Estimated $\hat{m}$ (stationary MR)")
    ax.set_title("Stationary MR is non-injective\nnear criticality", fontsize=8.5)
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    save(fig, "Figure3_MR_fold")


def fig4():
    d = pd.read_csv(f"{PROJ}/transient_onset_validation.csv")
    full = d[d.subsample_f == 1.0].sort_values("m_true")
    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    ax.plot([0.75, 1.35], [0.75, 1.35], color=C["grey"], lw=0.8, ls=":", label="identity")
    ax.axhline(1.0, color=C["grey"], lw=0.7, ls=":")
    ax.axvline(1.0, color=C["grey"], lw=0.7, ls=":")
    ax.errorbar(full.m_true, full.m_onset, yerr=full.m_onset_sd, fmt="o-",
                color=C["green"], lw=1.4, ms=5, capsize=2, label="transient-onset")
    ax.plot(full.m_true, full.m_MR_stationary, "s--", color=C["blue"], lw=1.2, ms=4,
            label="stationary MR (same traces)")
    ax.annotate("crosses 1\nat $m$=1.00", xy=(1.0, 1.011), xytext=(0.82, 1.08), fontsize=7,
                arrowprops=dict(arrowstyle="->", color=C["green"], lw=0.9), color=C["green"])
    ax.set_xlabel("True generative $m$"); ax.set_ylabel(r"Estimated $\hat{m}$")
    ax.set_title("Transient-onset estimator recovers\nthe sign of $(m-1)$", fontsize=8.5)
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    save(fig, "Figure4_transient_onset")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
    print("all figures written")
