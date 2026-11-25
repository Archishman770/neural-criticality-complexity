TITLE: How Bin-Width Choices Manufacture Structure in Branching-Ratio Estimates of Neural Criticality
AUTHOR: Archishman Sai Vemula
AFFILIATION: Texas Academy of Mathematics and Science, University of North Texas, Denton, TX 76203, USA
CORRESPONDENCE: Correspondence: [AUTHOR: corresponding e-mail]

ABSTRACT: Statistical complexity is hypothesised to peak near neural criticality, where the branching ratio m equals 1; testing this requires reliable estimates of m. Two routine, bin-specified analysis choices can manufacture structure in branching-ratio analyses while passing standard robustness checks. In mouse Neuropixels spiking, a multistep-regression fit window fixed in bins imposed a ceiling on m that was identical at every bin width (maximum 0.979 = exp(−1/48)) and censored 28% of region-level estimates; a window fixed in physical time removed the censoring and raised the maximum to 0.993. Synthetic autoregressive processes with known timescales reproduced the censoring boundary. In propofol-sedation EEG, setting the bin width from the event rate placed event rate on both sides of a complexity–criticality regression; with the bin width fixed, the statistical-complexity effect weakened by about three orders of magnitude in p and did not survive an amplitude control. We further show that stationary multistep regression cannot return supercritical values on stationary data, and that a transient-onset estimator and the avalanche measures κ and DCC did not locate m = 1 under realistic subsampling in a bursting network model. Fit windows and bin widths should be fixed in physical time, and criticality estimators validated with observed-only protocols.

LEAD: [AUTHOR: Chaos lead paragraph; delete if the submission system does not request one] Many complex systems are thought to operate near a critical point, where activity neither dies out nor runs away. In neural recordings this is usually judged by estimating the branching ratio, the average number of events triggered by each event, which equals one at criticality. We show that two ordinary choices made while estimating it, the length of a fit window counted in time bins and a bin width set from the event rate, can create apparent results that standard robustness checks confirm instead of catching. We also show that, under the partial sampling typical of real recordings, none of the estimators we tested located the critical point. Fixing analysis windows and bin widths in physical time removes the two artefacts.

KEYWORDS: branching ratio; criticality; multistep regression; statistical complexity; Lempel–Ziv complexity; neuronal avalanches; subsampling; estimator bias; EEG; Neuropixels

# I. INTRODUCTION

The hypothesis that cortical networks operate near a critical point has shaped two decades of work on neural dynamics. Neuronal avalanches with power-law size distributions were its first widely cited signature [@beggs2003], and criticality has since been linked to maximal dynamic range [@shew2009]. A more direct order parameter is the branching ratio m, the expected number of events triggered by one event: activity damps for m < 1, amplifies for m > 1, and m = 1 is the critical point. The multistep-regression (MR) estimator infers m from the decay of the population autocorrelation and is robust to the spatial subsampling that biases avalanche statistics [@wilting2018; @priesemann2009]. Applied to spiking in several species, it places cortex in a slightly subcritical, "reverberating" regime [@priesemann2014; @wilting2018; @wilting2019].

A parallel literature relates the complexity of neural signals to criticality. Statistical complexity in the sense of Rosso and colleagues, computed from ordinal-pattern distributions [@bandt2002; @rosso2007; @lamberti2004], vanishes for both perfectly ordered and perfectly random signals and is maximal in between, whereas Lempel–Ziv (LZ) complexity [@lempel1976] increases monotonically with randomness. Neural complexity peaks near criticality in cultured networks [@timme2016], statistical complexity is maximal close to a cortical phase transition in anaesthetised rats [@lotfi2021], and in human EEG, distance from criticality and complexity-related measures together track anaesthetic-induced loss of consciousness [@maschke2024; @casali2013; @schartner2015]. If complexity peaks near criticality, the location of that peak relative to m = 1, and whether it is the same across preparations, becomes a natural quantitative question. Answering it requires an estimate of m that is trustworthy on both sides of the critical point.

We set out to measure that offset in openly available human EEG and mouse spiking data and found instead that the estimate of m depended on analysis choices in ways the field's standard checks do not reveal. This paper reports those findings. It makes four contributions:

- A fit-quality window specified in bins imposes, through the identity m = exp(−Δt/τ), a ceiling on m that is identical at every bin width, so a bin-width robustness check cannot detect it (Sec. III A).
- Setting the bin width from the event rate places event rate on both sides of a complexity–criticality regression, a confound that within-subject label permutation cannot detect (Sec. III B).
- Stationary MR cannot return supercritical values on stationary data, and in a saturating branching model its estimate folds back below 1 for supercritical generative m (Sec. III D).
- A transient-onset estimator and the avalanche measures κ and DCC, tested under observation protocols available in real recordings, did not locate m = 1 in a subsampled bursting network (Sec. III E).

The neuroscience in the EEG and spiking arms is used as worked examples, not as the claims of the paper.

# II. MATERIALS AND METHODS

## A. Data

Human EEG came from a propofol-sedation study [@chennu2016; @chennudata]: 20 healthy adults recorded at baseline, mild sedation, moderate sedation and recovery, about 7 min of eyes-closed resting EEG per level. We used the data as distributed, with 91 scalp channels at 250 Hz, filtered at 0.5–45 Hz, cleaned, average-referenced and segmented into 10 s epochs by the original authors, giving 80 subject-condition cells. The four levels were recorded in a fixed sequence.

Mouse spiking came from the Allen Brain Observatory Visual Coding Neuropixels dataset [@siegle2021]. Four sessions were analysed (three functional-connectivity sessions and one Brain Observatory 1.1 session), using only spontaneous grey-screen intervals. Units were retained if labelled "good" with amplitude cutoff ≤ 0.1, presence ratio ≥ 0.9 and ISI violations ≤ 0.5, and spikes were pooled within brain regions having at least 20 retained units, giving 53 session–region estimates.

## B. Complexity measures

LZ complexity was computed on each channel and epoch after binarising at the median, using LZ76 phrase counts normalised by n/log₂n [@lempel1976]. Statistical complexity, here MPR_C, was computed from Bandt–Pompe ordinal patterns of order d = 6 with delay 1 [@bandt2002], as the product of normalised permutation entropy and the Jensen–Shannon disequilibrium normalised as in [@lamberti2004; @rosso2007]. Both implementations were checked against synthetic signals before use: on a sine wave, Gaussian noise and the logistic map at r = 4 (n = 5000), LZ returned 0.012, 1.037 and 1.025, so noise and chaos are nearly indistinguishable to it, whereas MPR_C returned 0.211, 0.028 and 0.484. These values are pinned by regression tests in the code repository. Cell values were averaged over channels and epochs.

## C. Population activity and binning

For EEG, each channel was z-scored per epoch, events were the peaks of |z| exceeding 2.5 SD, and event counts were summed across channels. In the original pipeline the summed series was binned at the mean inter-event interval of each cell, clamped to 1–12 samples. The realised range was 5–11 samples (20–44 ms), so the bin width Δt is approximately the reciprocal of the event rate. The binning-controlled analysis re-estimated m for every cell at each fixed Δt from 5 to 11 samples. For spiking data, population activity was the summed spike count per region in fixed bins of 4, 8 or 16 ms.

## D. Multistep-regression estimation and the fit-quality window

The branching ratio was estimated with the MR estimator (MR. Estimator toolbox, trial-separated coefficients, exponential fit) [@wilting2018; @spitzner2021]. The fit returns the autocorrelation time τ and m = exp(−Δt/τ), so m and τ are one quantity reparameterised and m depends on Δt. A fit was accepted if it converged, if 2 bins ≤ τ ≤ 0.8·KMAX, where KMAX is the largest lag fitted, and if 0 < m < 1.2 (EEG) or 0 < m < 1.5 (spiking). Through the identity above, the constraint τ ≤ 0.8·KMAX is also a constraint m ≤ exp(−1/(0.8·KMAX)). The original pipeline used KMAX = 60 bins in both arms. The time-specified alternative fixes KMAX in milliseconds: for spiking we report KMAX from 60 to 960 bins, and the window of 3840 ms (960, 480 and 240 bins at 4, 8 and 16 ms). For EEG the 60-bin window corresponds to a median of 1680 ms; we also swept it from 960 to 7680 ms. Distance from criticality was dist = |1 − m|.

## E. Statistics

Complexity was regressed on dist with linear mixed-effects models with a per-subject random intercept, measure ~ dist + (1 | subject). Every such model was refitted with raw signal variance as a covariate, as a control for amplitude. Significance of the complexity–dist relationship was also assessed with 10,000 within-subject permutations of condition labels, using the within-subject correlation between complexity and dist as the statistic. Comparisons between LZ and MPR_C used z-scored outcomes, because the two measures differ in scale by a factor of 3.5. Differences in m across the four conditions were tested with Friedman tests, and differences in LZ from baseline with Wilcoxon signed-rank tests.

## F. Synthetic models

Three generative models with known ground truth were used.

- **AR(1) processes** with autocorrelation time τ ∈ {40, 100, 200} steps (12 trials of 54,000 steps), summed into bins of 1, 2 or 4 steps. For these, the MR identity m = exp(−1/τ) holds exactly.
- **A saturating branching process**: A(t + 1) ~ Poisson(m·A(t)·(1 − A(t)/N) + h), with N = 1000, h = 2, m from 0.8 to 1.3, and observation fractions f = 1 and 0.1 (binomial thinning).
- **A bursting network with synaptic depression**: A(t + 1) ~ Poisson(m·r(t)·A(t)·(1 − A(t)/N) + h), with resources r(t + 1) = r + (1 − r)/τ_rec − u·A(t + 1)/N, N = 1000, and seven parameter settings (h ∈ {0.1, 0.2, 0.4}, τ_rec ∈ {150, 300, 600} steps, u ∈ {0.03, 0.05, 0.08}, varied one at a time). Each setting was run for 100,000 steps per seed with 4 seeds, and observed with binomial thinning at f ∈ {1, 0.10, 0.05, 0.02}.

## G. Estimators under subsampling

Three estimators were evaluated on the bursting network under protocols that use only what a recording provides.

- **Transient-onset estimator.** Burst onsets are aligned, and m is estimated from the growth of the onset-aligned ensemble mean, m̂ = exp(slope of log⟨A⟩ against time), over the unsaturated window. We compared three detection protocols:
  - onsets detected on the full simulated population ("oracle", as in the estimator's original validation);
  - onsets detected on the observed counts, with a threshold of 2 or 3 observed events after 6 quiet bins and a fit ceiling of 10% of the median observed burst peak;
  - a split-half protocol, detecting onsets on one random half of the observed units and measuring growth on the other.
- **The avalanche measure κ** [@shew2009], with a critical size exponent of 1.5 and 10 log-spaced comparison points.
- **The deviation from criticality coefficient (DCC)** [@ma2019; @fontenele2019], computed from size and duration exponents fitted by discrete power-law maximum likelihood with a Kolmogorov–Smirnov choice of the lower cutoff [@clauset2009].

Avalanches were runs of non-empty bins at 1, 2 and 4 model steps per bin.

Rank stability of κ was assessed on 400 synthetic recordings. Each had m ~ U(0.8, 1.3), network parameters drawn from the ranges above, and a fixed random set of 20–142 observed units with log-normal firing propensities, so that the effective observation fraction differed across recordings (0.008–0.092). Each recording was simulated for 300,000 steps after a 5,000-step burn-in, with κ computed at 2, 4 and 8 steps per bin.

## H. Software

Analyses used Python with NumPy, SciPy, pandas, statsmodels, MR. Estimator [@spitzner2021], MNE and pymatreader. Every numerical result in this paper is produced by a numbered script in the code repository, and every figure is generated from those scripts' outputs.

# III. RESULTS

## A. A fit window fixed in bins imposes a bin-invariant ceiling on m

With the original 60-bin window, the largest branching ratio in the spiking data was 0.979 at 4 ms bins, which equals the filter ceiling exp(−1/48) = 0.97938 to three decimals (Table 1, Figure 1A). The filter rejected 28.3% of session–region estimates at 4 ms as having no decay within the fit window, which is exactly where near-critical estimates would lie (Figure 1B). Because the ceiling is set in bins, it is the same at every bin width. A bin-width sweep, the usual robustness check, therefore returns an apparently stable ceiling, because the check and the artefact share the same invariance.

Lengthening the window removed the ceiling. At 4 ms, censoring fell from 28.3% at 60 bins to 13.2% at 120 bins and to zero from 240 bins upward. The maximum m rose to 0.993 and the median to 0.978, with 22 of the 47 passing estimates at or above 0.98 (none at 60 bins). The estimates had converged: median m changed by less than 0.001 between 240 and 960 bins. The fitted autocorrelation time at 4 ms rose from a median of 71 ms to 177 ms. At 16 ms bins the 60-bin window already spans 960 ms, and estimates were insensitive to the window (maximum 0.975–0.977, median 0.939, median τ 253–254 ms).

**Table 1.** Branching-ratio estimates for the spiking data (53 session–region estimates) under fit windows specified in bins and in time.

| Bin | Fit window | Passing | Censored | Max m | Median m | Median τ (ms) | Estimates ≥ 0.98 |
|---|---|---|---|---|---|---|---|
| 4 ms | 60 bins (240 ms) | 32 | 28.3% | 0.979 | 0.945 | 71 | 0 |
| 4 ms | 120 bins (480 ms) | 40 | 13.2% | 0.989 | 0.972 | 140 | 13 |
| 4 ms | 960 bins (3840 ms) | 47 | 0% | 0.993 | 0.978 | 177 | 22 |
| 16 ms | 60 bins (960 ms) | 42 | 0% | 0.977 | 0.939 | 253 | 0 |
| 16 ms | 240 bins (3840 ms) | 43 | 0% | 0.975 | 0.939 | 254 | 0 |

FIGURE: figures_entropy/Fig1_fit_window_ceiling.png | **Figure 1.** The fit-quality window sets a ceiling on m. (**A**) Largest branching ratio across 53 session–region estimates as a function of the fit window KMAX, for 4 ms and 16 ms bins. The dashed curve is the ceiling exp(−1/(0.8·KMAX)) implied by the constraint τ ≤ 0.8·KMAX. At 4 ms the maximum follows the ceiling until the window is long enough, then levels off near 0.993. (**B**) Percentage of estimates censored by the constraint. No 16 ms estimate passed at 960 bins, because recordings were shorter than the required trial length.

Synthetic AR(1) processes with known τ isolate the mechanism. Under the 60-bin rule, every process with true τ ≤ 40 bins passed. Four of the five processes with true τ ≥ 50 bins were censored; the fifth (true τ = 50, fitted τ = 46.7) passed narrowly, because censoring acts on the fitted rather than the true τ. Fitting with the longer window changed the recovered m by at most 0.001, and the time-specified estimates recovered the true m to within 0.007. On a single-timescale process, then, a short window censors but does not bias the magnitude. The roughly 2.5-fold increase in fitted τ in the spiking data when the window was lengthened therefore indicates autocorrelation with more than one timescale, which the short window under-fits.

The EEG arm was not affected in this way. Its 60-bin window spanned about 5.5 times the median τ (297 ms), no cell was censored, and sweeping the window from 960 to 7680 ms changed the median m by 0.0001 and the median τ by 3%.

## B. Setting the bin width from the event rate confounds a complexity–criticality regression

In the EEG data, statistical complexity appeared to increase strongly as m approached 1. With the original binning, MPR_C ~ dist + (1 | subject) gave β = −0.154 (p = 6.4 × 10⁻⁶), and no within-subject permutation of condition labels reproduced a correlation that large (0 of 10,000). This relationship has the binning parameter on both sides. Δt was set to the mean inter-event interval, so it is approximately the reciprocal of the event rate; it enters m through m = exp(−Δt/τ), while the event rate enters the complexity measures. Across cells, LZ correlated with events per epoch at r = +0.94 (Figure 2A) and with Δt at r = −0.80, and MPR_C correlated with Δt at r = −0.72. Within-subject permutation cannot detect this, because the confound is carried by both permuted variables rather than aligned with condition labels.

Holding Δt fixed removes the shared path (Table 2, Figure 2B,C). At each fixed Δt from 5 to 11 samples, the MPR_C effect remained negative but was about three orders of magnitude weaker (β = −0.151 to −0.075; p = 0.013–0.024; permutation p = 0.007–0.011). It did not survive the amplitude control: adding raw variance removed 66–71% of the coefficient and left p = 0.40–0.50 at every Δt. Under the original binning, the same control had appeared to pass (p = 0.010, 37% attenuation). LZ kept most of its coefficient with variance in the model (7–12% attenuation), remaining significant at four of seven bin widths (p = 0.016–0.077). Adjusting the original model for events per epoch, rather than fixing Δt, also removed the effect (LZ β = +0.075, p = 0.155; MPR_C β = −0.037, p = 0.244). We prefer the design-based control, because Δt alone accounts for about 64% of the LZ variance and regressing it out removes much of the outcome along with the confound. On a scale-free measure, the within-subject correlation between complexity and dist fell from |r| = 0.47–0.48 to 0.28–0.35.

**Table 2.** Complexity–criticality models in the EEG data (80 cells, 20 subjects) under the original and fixed binning. β and p are from the mixed model; the permutation p uses 10,000 within-subject permutations of the complexity–dist correlation. Ranges are across Δt = 5–11 samples.

| Binning | Measure | β (dist) | p | Permutation p | p with variance | Attenuation |
|---|---|---|---|---|---|---|
| Event-rate (original) | MPR_C | −0.154 | 6.4 × 10⁻⁶ | ≤ 10⁻⁴ | 0.010 | 37% |
| Event-rate (original) | LZ | −0.479 | 6.9 × 10⁻⁶ | 6 × 10⁻⁴ | 2.2 × 10⁻⁵ | −11% |
| Fixed Δt | MPR_C | −0.151 to −0.075 | 0.013–0.024 | 0.007–0.011 | 0.40–0.50 | 66–71% |
| Fixed Δt | LZ | −0.548 to −0.214 | 0.005–0.027 | 0.011–0.022 | 0.016–0.077 | 7–12% |

FIGURE: figures_entropy/Fig2_binning_confound.png | **Figure 2.** A bin width set from the event rate confounds the complexity–criticality relationship in propofol EEG. (**A**) LZ complexity against events per epoch across the 80 subject-condition cells. (**B**, **C**) −log₁₀ p of the dist coefficient for MPR_C (**B**) and LZ (**C**), with event-rate binning ("data-driven") and at each fixed bin width, without (circles) and with (triangles) raw variance as a covariate. The dashed line marks p = 0.05.

## C. What remains in the EEG once binning is controlled

With Δt fixed, dist is a monotone transform of τ (r = −0.89 to −0.91 across Δt), and complexity tracks τ directly (LZ p = 0.0008–0.0074; MPR_C p = 0.018–0.053). For a randomness measure such as LZ, lower complexity at longer autocorrelation time is close to definitional, so the surviving relationship is not evidence specific to criticality.

The branching ratio itself did not differ across the four conditions (Friedman χ² = 1.50, p = 0.68 with event-rate binning; p = 0.66–0.88 at each fixed Δt). LZ was higher than at baseline in all three later conditions (mild +0.025, 17 of 20 subjects, p = 0.005; moderate +0.027, p = 0.040; recovery +0.022, p = 0.003). An increase at low propofol doses has been reported on this dataset [@boncompte2021], in contrast to the decrease under general anaesthesia [@schartner2015]. That LZ remained elevated at recovery, however, points to an order or baseline-recording contribution, which the fixed recording sequence does not allow to be separated from a drug effect. These sedation data do not reach the anaesthetic-induced loss of consciousness at which a shift in branching ratio has been reported [@maschke2024]. We therefore claim neither a replication nor a contradiction of that finding.

## D. Stationary MR cannot return supercritical values on stationary data

The MR regression slope at lag k is Cov(A_t, A_{t+k})/Var(A_t). For a stationary process this is the lag-k autocorrelation, which cannot grow without bound, so the fitted exponential decays (τ > 0) and m = exp(−Δt/τ) < 1, up to estimation noise. The bound comes from stationarity, not from the exponential form. It is implicit in the driven branching model on which MR rests, where a stationary mean activity h/(1 − m) exists only for m < 1 [@wilting2018].

The saturating branching model shows what this means for supercritical systems (Figure 3). Stationary MR rose with the generative m up to m̂ = 0.909 at m = 1.0 and then folded back, to 0.896 at m = 1.05 and 0.679 at m = 1.3. A strongly supercritical system thus reads as more subcritical than a subcritical one, and m̂ ≈ 0.89 is returned for both m = 0.95 and m = 1.05. The estimator follows the local slope of the dynamics at the stable operating point, which is m below criticality and 2 − m above it, rather than the generative m. Observing 10% of units changed the estimates by less than 0.02, so the fold is not a subsampling effect.

FIGURE: figures_entropy/Fig3_MR_fold.png | **Figure 3.** Stationary MR is non-injective around criticality. MR estimates against the generative branching ratio of a saturating branching process (N = 1000) with full (circles) and 10% (squares) observation. The dotted line is the identity. The dashed line is the local slope of the dynamics at the stable fixed point, which MR follows instead of the generative m.

## E. Estimators that target the supercritical side, under realistic subsampling

Because stationary MR cannot see supercriticality, we tested three alternatives on the bursting network, using only information available in a recording (Figure 4).

**Transient-onset estimator.** With onsets detected on the full simulated population, as in its original validation, the estimate crossed 1 near the true critical point (crossing at m = 0.94–1.07 across the seven parameter settings and three observation fractions). With onsets detected on the observed counts at threshold 2, the crossing moved to m = 1.04–1.28, and in 5 of the 21 conditions there was no crossing at all. The sign of (m − 1) was then correct for 38–100% of the generative values, depending on the setting, and for 38–71% at 2% observation. With threshold 3, three or four of the eight values had no estimate at 2% observation, and at most one of the rest was classified correctly. Split-half detection did not remove the bias (crossing 0.89–1.28). In the base setting at m = 1, the median true population activity at observed onsets was only 7–13 units, depending on the observation fraction (versus 3 for the oracle threshold), well below saturation. This suggests that many observed "onsets" are fluctuations that do not grow, although we have not shown this to be the only cause.

**κ.** κ increased with m in 80 of 84 combinations of parameter setting, observation fraction and bin width; all four exceptions were at 2% observation. Its crossing of 1 did not mark m = 1: it ranged from 0.82 to 1.29, or did not occur at all, depending on the observation fraction and bin width. With full observation κ exceeded 1 even at m = 0.8 (κ = 1.012), and at 2% observation it stayed below 1 even at m = 1.3 (κ = 0.888) (Figure 4B). The location of the DCC minimum ranged over the whole tested interval (m = 0.80–1.30), so DCC carried no information about m in this model.

**Ranking across heterogeneous recordings.** Across 400 recordings with different network parameters and observed units, κ at 4 steps per bin ranked m well overall (Spearman ρ = 0.86; 0.89–0.94 within each tercile of observation fraction; 0.64 for |m − 1| ≤ 0.1; 6.7% of pairs separated by Δm ≥ 0.1 ranked in reverse). κ was, however, almost a monotone function of the observed firing rate (Spearman ρ = 0.991; Figure 4C). After flexible adjustment for rate (a natural cubic spline with four degrees of freedom), 2.5% of the variance of κ remained, and κ added 0.021 to the variance in m explained by rate. In this model, κ orders m largely because both track firing rate.

FIGURE: figures_entropy/Fig4_subsampled_estimators.png | **Figure 4.** Criticality estimators under realistic subsampling in a bursting network with synaptic depression. (**A**) Transient-onset estimates at 5% observation with onsets detected on the full population (circles), on the observed units (triangles), or on one half of the observed units with growth measured on the other (squares). (**B**) κ against the generative m at four observation fractions f. (**C**) κ against the observed firing rate across 400 synthetic recordings with heterogeneous parameters and observed units.

**Table 3.** Summary of estimators of criticality tested under observation protocols available in real recordings (bursting network, seven parameter settings).

| Estimator | Reads m > 1 | Locates m = 1 under realistic subsampling | Main failure |
|---|---|---|---|
| Stationary MR | No | Subcritical side only | Folds back for m > 1 (Sec. III D) |
| Transient onset, oracle onsets | Yes | Crossing 0.94–1.07 | Requires unobservable onsets |
| Transient onset, observed onsets | Partly | Crossing 1.04–1.28 or absent | Biased toward subcritical readings |
| κ | Yes | κ = 1 at m = 0.82–1.29 or never | Calibration shifts with observation fraction and bin width; tracks firing rate |
| DCC | — | Minimum at m = 0.80–1.30 | Uninformative in this model |

# IV. DISCUSSION

Two routine choices in branching-ratio analysis can manufacture structure that the usual checks do not detect. A fit-quality window specified in bins turns a bound on τ into a bound on m that is identical at every bin width, so a bin-width sweep confirms the artefact instead of exposing it. A bin width set from the event rate puts the event rate on both sides of any regression of a rate-sensitive quantity on m, and within-subject label permutation cannot reveal this because the confound moves with both variables. In both cases the robustness check shares the invariance of the artefact. The recommendations follow directly:

- specify MR fit windows in physical time and confirm that the fitted τ lies well inside them;
- hold the bin width fixed in physical time across conditions whenever m is related to other signal statistics;
- report τ, or m together with its bin width, since m alone is not a property of the tissue.

The limits we report on estimating supercritical m are related but distinct. Stationary MR is bounded below 1 on stationary data, so the absence of supercritical MR estimates in stationary recordings is a property of the method applied to such data, not a finding about cortex. This is consistent with the reverberating-regime picture [@wilting2018; @wilting2019], and our contribution here is to state it as an identifiability limit and to show the fold in a saturating model. The alternative estimators we tested did not escape it once evaluated with the information a recording provides. The transient-onset estimator worked only when onsets were taken from the unobserved full population. κ ranked m within a fixed observation fraction but did not locate m = 1, which is consistent with earlier demonstrations that subsampling distorts avalanche statistics [@priesemann2009]. In our heterogeneous population, κ also tracked firing rate closely. We stress that these are results for three estimator implementations on one model family with binomial or propensity-weighted thinning. They do not show that absolute criticality is unidentifiable in real data, and estimators designed for spatially structured subsampling may behave differently.

For the question that motivated this work, the location of the complexity peak relative to m = 1, the implications are practical. The peak cannot be located on the m axis with an estimator that cannot read m > 1, or whose crossing of 1 moves with sampling. Previous reports that complexity is maximal near criticality [@timme2016; @lotfi2021] used other order parameters or relative measures of criticality, and our results do not bear on their conclusions. They indicate, however, that a cross-preparation comparison of where complexity peaks on the m axis requires an estimator validated to locate m = 1 under the sampling conditions of each preparation.

The EEG example also carries a caution for complexity–criticality analyses in general. Once binning was controlled, the relationship in these data reduced largely to complexity covarying with autocorrelation time, and the statistical-complexity component could not be separated from slow-wave amplitude. The branching ratio did not change across sedation levels.

This study has limitations. The EEG population signal was built from thresholded scalp events, and branching ratios derived this way are not equivalent to spiking branching ratios. The spiking arm comprised four sessions of mixed type. The synthetic models omit spatial structure, neuron-level refractoriness and heterogeneous burst participation. The parameters of the transient-onset and avalanche analyses (onset thresholds, quiet windows, bin widths) were examined over limited ranges.

# V. CONCLUSIONS

Bin-specified choices in branching-ratio estimation can create a subcritical ceiling and a complexity–criticality relationship that standard robustness checks confirm rather than expose. Fixing fit windows and bin widths in physical time removes both. On stationary data, multistep regression cannot report supercritical values, and the transient-onset and avalanche alternatives we tested did not locate the critical point under realistic subsampling. Questions that depend on where a quantity peaks relative to m = 1 therefore need estimators validated with observed-only protocols before they can be answered.

BACKMATTER:

# ACKNOWLEDGMENTS

[AUTHOR: confirm funding] This research received no external funding. [AUTHOR: confirm wording] During the preparation of this work the author used Claude (Anthropic) to assist with code development, analysis and drafting. The author reviewed and edited all content and takes full responsibility for the publication.

# AUTHOR DECLARATIONS

**Conflict of Interest:** The author has no conflicts to disclose.

**Author Contributions:** **Archishman Sai Vemula:** Conceptualization (lead); Methodology (lead); Software (lead); Formal analysis (lead); Investigation (lead); Visualization (lead); Writing – original draft (lead); Writing – review & editing (lead).

**Ethics Approval:** This study is a secondary analysis of openly available, de-identified data. The human EEG study was approved by the Cambridgeshire 2 Regional Ethics Committee and conducted in accordance with the Declaration of Helsinki, with written informed consent from all participants [@chennu2016]. Animal procedures for the mouse recordings were approved by the Allen Institute's Institutional Animal Care and Use Committee [@siegle2021]. No new data were collected.

# DATA AVAILABILITY

The EEG data are openly available in the University of Cambridge Repository at https://doi.org/10.17863/CAM.68959 [@chennudata]. The mouse Neuropixels data are openly available from the Allen Brain Observatory at https://portal.brain-map.org/explore/circuits/visual-coding-neuropixels [@siegle2021]. The code that supports the findings of this study, including all synthetic-data generators and figure scripts, is openly available at [AUTHOR: insert public repository URL and archived DOI].

REFERENCES
beggs2003 | J. M. Beggs and D. Plenz, "Neuronal avalanches in neocortical circuits," J. Neurosci. **23**, 11167–11177 (2003).
shew2009 | W. L. Shew, H. Yang, T. Petermann, R. Roy, and D. Plenz, "Neuronal avalanches imply maximum dynamic range in cortical networks at criticality," J. Neurosci. **29**, 15595–15600 (2009).
wilting2018 | J. Wilting and V. Priesemann, "Inferring collective dynamical states from widely unobserved systems," Nat. Commun. **9**, 2325 (2018).
priesemann2009 | V. Priesemann, M. H. J. Munk, and M. Wibral, "Subsampling effects in neuronal avalanche distributions recorded in vivo," BMC Neurosci. **10**, 40 (2009).
priesemann2014 | V. Priesemann, M. Wibral, M. Valderrama, R. Pröpper, M. Le Van Quyen, T. Geisel, J. Triesch, D. Nikolić, and M. H. J. Munk, "Spike avalanches in vivo suggest a driven, slightly subcritical brain state," Front. Syst. Neurosci. **8**, 108 (2014).
wilting2019 | J. Wilting and V. Priesemann, "Between perfectly critical and fully irregular: A reverberating model captures and predicts cortical spike propagation," Cereb. Cortex **29**, 2759–2770 (2019).
bandt2002 | C. Bandt and B. Pompe, "Permutation entropy: A natural complexity measure for time series," Phys. Rev. Lett. **88**, 174102 (2002).
rosso2007 | O. A. Rosso, H. A. Larrondo, M. T. Martín, A. Plastino, and M. A. Fuentes, "Distinguishing noise from chaos," Phys. Rev. Lett. **99**, 154102 (2007).
lamberti2004 | P. W. Lamberti, M. T. Martin, A. Plastino, and O. A. Rosso, "Intensive entropic non-triviality measure," Physica A **334**, 119–131 (2004).
lempel1976 | A. Lempel and J. Ziv, "On the complexity of finite sequences," IEEE Trans. Inf. Theory **22**, 75–81 (1976).
timme2016 | N. M. Timme, N. J. Marshall, N. Bennett, M. Ripp, E. Lautzenhiser, and J. M. Beggs, "Criticality maximizes complexity in neural tissue," Front. Physiol. **7**, 425 (2016).
lotfi2021 | N. Lotfi, T. Feliciano, L. A. A. Aguiar, T. P. L. Silva, T. T. A. Carvalho, O. A. Rosso, M. Copelli, F. S. Matias, and P. V. Carelli, "Statistical complexity is maximized close to criticality in cortical dynamics," Phys. Rev. E **103**, 012415 (2021).
maschke2024 | C. Maschke, J. O'Byrne, M. A. Colombo, M. Boly, O. Gosseries, S. Laureys, M. Rosanova, K. Jerbi, and S. Blain-Moraes, "Critical dynamics in spontaneous EEG predict anesthetic-induced loss of consciousness and perturbational complexity," Commun. Biol. **7**, 946 (2024).
casali2013 | A. G. Casali, O. Gosseries, M. Rosanova, M. Boly, S. Sarasso, K. R. Casali, S. Casarotto, M.-A. Bruno, S. Laureys, G. Tononi, and M. Massimini, "A theoretically based index of consciousness independent of sensory processing and behavior," Sci. Transl. Med. **5**, 198ra105 (2013).
schartner2015 | M. Schartner, A. Seth, Q. Noirhomme, M. Boly, M.-A. Bruno, S. Laureys, and A. Barrett, "Complexity of multi-dimensional spontaneous EEG decreases during propofol induced general anaesthesia," PLoS ONE **10**, e0133532 (2015).
chennu2016 | S. Chennu, S. O'Connor, R. Adapa, D. K. Menon, and T. A. Bekinschtein, "Brain connectivity dissociates responsiveness from drug exposure during propofol-induced transitions of consciousness," PLoS Comput. Biol. **12**, e1004669 (2016).
chennudata | S. Chennu, S. O'Connor, R. Adapa, D. K. Menon, and T. A. Bekinschtein, "Research data supporting 'Brain connectivity during propofol sedation'," University of Cambridge Repository (2016), https://doi.org/10.17863/CAM.68959.
siegle2021 | J. H. Siegle, X. Jia, S. Durand, S. Gale, C. Bennett, N. Graddis, G. Heller, T. K. Ramirez, H. Choi, J. A. Luviano *et al.*, "Survey of spiking in the mouse visual system reveals functional hierarchy," Nature **592**, 86–92 (2021).
spitzner2021 | F. P. Spitzner, J. Dehning, J. Wilting, A. Hagemann, J. P. Neto, J. Zierenberg, and V. Priesemann, "MR. Estimator, a toolbox to determine intrinsic timescales from subsampled spiking activity," PLoS ONE **16**, e0249447 (2021).
ma2019 | Z. Ma, G. G. Turrigiano, R. Wessel, and K. B. Hengen, "Cortical circuit dynamics are homeostatically tuned to criticality in vivo," Neuron **104**, 655–664.e4 (2019).
fontenele2019 | A. J. Fontenele, N. A. P. de Vasconcelos, T. Feliciano, L. A. A. Aguiar, C. Soares-Cunha, B. Coimbra, L. Dalla Porta, S. Ribeiro, A. J. Rodrigues, N. Sousa *et al.*, "Criticality between cortical states," Phys. Rev. Lett. **122**, 208101 (2019).
clauset2009 | A. Clauset, C. R. Shalizi, and M. E. J. Newman, "Power-law distributions in empirical data," SIAM Rev. **51**, 661–703 (2009).
boncompte2021 | G. Boncompte, V. Medel, L. I. Cortínez, and T. Ossandón, "Brain activity complexity has a nonlinear relation to the level of propofol sedation," Br. J. Anaesth. **127**, 254–263 (2021).
