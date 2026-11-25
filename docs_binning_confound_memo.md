# Memo: the §3.1 binning confound and what survives it

Date: 2026-08-11
Scripts: `step19_fixed_dt_branching.py`, `step20_fixed_dt_test.py`
Outputs: `Sedation-RestingState/branching_summary_fixeddt.csv`,
`Sedation-RestingState/fixed_dt_test_results.csv`

## The problem

`step3b_branching_ratio_v2.py` line 42 sets the bin width from the data:

    dt = int(round(np.mean(np.concatenate(ieis))))

so dt is the mean inter-event interval, i.e. approximately 1/event-rate. Since
`dist = 1 − exp(−dt/τ)`, dt enters the criticality side directly. It also enters
the complexity side: in these data LZ correlates with `events_per_epoch` at
**r = +0.938** and with dt at **r = −0.801**; MPR_C with dt at **r = −0.723**.
The §3.1 relationship therefore has the binning parameter on both sides of it.

This was not caught by the existing robustness checks. Within-subject label
permutation (10,000 shuffles) cannot detect it, because the confound is shared by
both permuted variables rather than aligned with condition labels.

**Note:** `step4_join_and_test.py` already contains the relevant control as its
model 2. Re-running it gives LZ β_dist = +0.0752, **p = 0.156**; MPR_C
β_dist = −0.0373, **p = 0.244**, with events_per_epoch at p = 2×10⁻⁸⁸. That model
was computed but its result is not reported anywhere in the manuscript.

## The test

Held dt constant across all 80 subject-condition cells and re-estimated the
branching ratio, sweeping dt = 5…11 samples (the range step3b selects). Everything
else identical: 2.5-SD peak detection, trial-separated MR coefficients, 60-bin fit
window, same quality criteria. All 560 cells converged.

A design-based control (hold dt fixed) is preferable here to a statistical one
(regress dt out), because dt explains ~64% of LZ variance and adjusting for it
removes most of the outcome's variance along with the confound.

## What survives

| binning | measure | β_dist | p | p_perm (10k) | p after variance | attenuation |
|---|---|---|---|---|---|---|
| data-driven (manuscript) | LZ | −0.479 | 7×10⁻⁶ | 0.0006 | 2×10⁻⁵ | −11% |
| data-driven (manuscript) | MPR_C | −0.154 | 6×10⁻⁶ | 0.0001 | 0.0102 | 37% |
| fixed dt = 5…11 | LZ | −0.55 … −0.21 | 0.005–0.027 | 0.011–0.022 | 0.016–0.077 | 7–12% |
| fixed dt = 5…11 | MPR_C | −0.15 … −0.075 | 0.013–0.024 | 0.007–0.011 | **0.40–0.50** | **66–71%** |

Standardised (within-subject |r|), the scale-free comparison: **0.47–0.48 under
data-driven binning, 0.28–0.35 under fixed binning.**

Under fixed binning the effect is present but roughly three orders of magnitude
weaker in p than reported. It tracks τ directly (LZ p = 0.0008–0.0074; MPR_C
p = 0.018–0.053 across dt).

## Consequences for the manuscript

1. **§3.1's headline numbers are inflated and must be restated.** β = −0.1544,
   p = 6.4×10⁻⁶ is not the effect after binning control. The fixed-dt numbers
   (p ≈ 0.01–0.02) are the defensible ones.
2. **Falsification criterion 2 FAILS for MPR_C.** At every fixed dt the variance
   control removes 66–71% of the coefficient and leaves p ≈ 0.40–0.50. The
   manuscript currently reports criterion 2 as passing (p = 0.0102). Hard rules 1
   and 5 require this be stated as a failure.
3. **The §3.2 asymmetry survives and strengthens.** LZ retains a criticality
   relation independent of amplitude (attenuation 7–12%); MPR_C does not
   (66–71%). The manuscript's claim that MPR_C carries an extra amplitude
   component is more strongly supported under proper binning control, not less.
4. **What survives is "complexity tracks autocorrelation time."** Since
   m = exp(−Δt/τ), with dt fixed, dist is τ rescaled (r = −0.89 to −0.91). The
   surviving claim reduces to complexity covarying with τ, which for a
   randomness measure is close to definitional and is not a criticality result.

## Separate findings from the same pass

- **The branching ratio does not track propofol.** Condition explains **0.4%** of
  m variance; Friedman χ² = 1.50, **p = 0.68** across baseline/mild/moderate/
  recovery (means 0.9034, 0.9073, 0.9026, 0.9009). Subject explains 32%.
- **LZ moves the wrong way under sedation.** Mild vs baseline Δ = +0.0245,
  p = 0.0049, 17/20 subjects increased; recovery remains elevated (p = 0.0027).
  Literature (Schartner et al. 2015) reports LZ decreasing under propofol. The
  failure to return at recovery points to a session-order effect, with baseline as
  the outlier condition, rather than pharmacology. This warrants a look at the
  baseline recordings before §3.6 is relied on.

## Recommended next step

The binning confound is a second, independent instance of the failure mode already
documented in §3.3: a bin-specified choice manufacturing structure that the
field's standard robustness checks (bin-width invariance, label permutation)
cannot see. Reported as a pair, that is a stronger and more general methods
result than §3.3 alone.
