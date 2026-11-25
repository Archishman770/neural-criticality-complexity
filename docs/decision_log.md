# Decision Log

Running log of autonomous default choices (per `CLAUDE.md` Section 4) and
pre-specified analysis decisions. Newest at bottom. Each entry: what was decided,
why, and whether it's routine or a flagged (hard-rule / plan / dataset) decision.

Dates are approximate to the working session; this log was reconstructed at the end
of the propofol+Allen phase, so earlier routine choices may be summarized.

---

**[routine] Branching estimator: v2 only.** `branching_summary_v2.csv` (Wilting-
Priesemann MR) is the sole valid branching file. v1 (naive avalanche sigma) is
irreparably subsampling-biased; retained only as an explicitly labeled biased
reference, never used in analysis.

**[routine] MR quality filter (fixed pre-inspection).** `mr_ok` = fit converges AND
2.0 ≤ tau_bins ≤ 0.8·KMAX AND 0 < m_MR < 1.2. Cells failing get a `reason` and are
dropped from correlational analyses but reported (drop rate + reasons).

---

## Propofol permutation phase

**[flagged: hard rule 3] step6 within-subject permutation test.** Statistic = mean
over subjects of within-subject Spearman(m_MR, measure); 10,000 within-subject
label shuffles; variance travels alongside (rule 2). Chosen because it's fast, runs
off the existing CSVs, and stress-tests the within-subject effect (the estimand most
exposed to small-n optimism). Explicitly logged as a **robustness check, not
confirmatory** — pipeline was post-hoc tuned. Results: LZ ρ=+0.390 p=0.0011; MPR_C
ρ=+0.440 p=0.0008; variance ρ=−0.260 p=0.0498.

**[note] Falsification pass/fail deferred.** step6 tests the within-subject Spearman,
not criterion 1's mixed-model `dist` coefficient. Explicit pass/fail on criteria 1–3
is owed against the step4/step5 mixed-model output and remains outstanding.

---

## Allen mouse Neuropixels phase (cross-species / cross-modality)

**[flagged: new dataset] Allen Visual Coding Neuropixels adopted.** Adds species +
modality axes. Chosen over sleep EEG / psychedelic EEG as the first extension
because dense spiking was the best chance to escape the EEG subsampling ceiling and
reach m≈1.

**[flagged: environment] Read NWB directly via h5py; abandon AllenSDK reader.**
AllenSDK 2.16.2 + pynwb/hdmf/ndx-events is version-broken on this machine and could
not be pinned into a working state. The NWB files are valid HDF5; h5py reads them
with zero AllenSDK in the path. Standing default for all spiking data.

**[routine] MR on spiking = population spike-count per bin.** a(t) = total
population spike count per BIN_MS bin per region; spontaneous windows as trials.
Same MR estimator + quality criteria as step3b.

**[routine → flagged] Running/stationary state split DROPPED for the diagnostic.**
Initially split spontaneous windows by running speed (RUN≥3, STAT≤0.5 cm/s). Awake
mice are rarely still → "stationary" was always empty, "running" fragmented sparse
data into un-fittable pieces (non-finite MR coefficients). Decision: **pool all
spontaneous windows per region** for the go/no-go m-scan. The arousal split is
**deferred to analysis stage**, not deleted; if used it gets its own sensitivity
check + a log entry. (This touches how a criticality estimate is computed, hence
flagged.)

**[routine] Unit QC for Allen.** quality=='good' AND amplitude_cutoff≤0.1 AND
presence_ratio≥0.9 AND isi_violations≤0.5; regions kept only with ≥20 good units.

**[flagged: hard rule / sensitivity] Bin-width sweep 4/8/16 ms.** The 4 ms default
materially affects m, so before finalizing the subcritical NO-GO we swept bin width.
Result: max m ≤ 0.977 at all three widths; 0 cells ≥0.98; 0 supercritical. Ceiling
is bin-robust → NO-GO on awake in-vivo is solid, not a binning artifact. This is the
sensitivity check that had been flagged as pending. Explicitly a robustness check,
exploratory, n=4 sessions.

**[note] Scan-only m ceiling widened to 1.5.** For the Allen scan's quality filter,
the upper m bound was raised from step3b's 1.2 to 1.5, because supercritical values
are the *signal* in this scan, not noise to reject. Scan-specific; NOT propagated to
the EEG pipeline.

---

## Behavioral split phase (propofol)

**[flagged: plan amendment — see research_design A4] Behavioral responsiveness
split.** Moderate condition only. Threshold **unresponsive = correct-responses ≤ 20**
(max 40, ≤50%), locked from the behavioral distribution BEFORE complexity was
inspected → 6 unresponsive / 14 responsive. Directional hypothesis unresponsive <
responsive for MPR_C, LZ, m_MR.

**[flagged: analysis choice] Concentration included as covariate.** Because conc at
moderate is genuinely titrated (range 433–1521), responsive-vs-unresponsive is not
just low-vs-high dose. Including conc lets the claim be "responsiveness predicts
complexity controlling for dose" — the consciousness-not-pharmacology framing. If
conc had been constant the covariate would be dropped (collinear); it wasn't.

**[result] SUPPORTED (exploratory).** Group permutation: MPR_C p_dir=0.0034, LZ
0.0014, m_MR 0.0049; variance null (p=0.073). Dose-adjusted: MPR_C p=0.021, LZ
0.0006, m_MR 0.045. Caveats logged: n=6 → exploratory only; m_MR borderline; LZ
strongest. Sensitivity variant ≤21 (7v13) available as behavior-only robustness.

---

## Outstanding / owed

- Explicit written falsification pass/fail (criteria 1–3) against the step4/step5
  mixed-model output for the propofol arm.
- Rename the on-disk project folder (drop "thermodynamic invariant").
- Keep `CLAUDE.md`, this log, and current output CSVs
  (`step6_permutation_results.csv`, `allen_m_scan_binsweep.csv`,
  `step8_behavioral_split.csv`) in the Project knowledge panel, not just the repo.

---

## Log gap (noted 2026-09-26)

This file was copied into `docs/` from `files/decision_log.md` (last edited
2026-07-05) to become the canonical log that CLAUDE.md Section 6 names. Amendments
A5–A7 and the steps 14–17 decisions are cited in the manuscript and scripts but were
never entered here; they are still owed (see the "Outstanding" list above).

---

## Binning confound phase (propofol) — amendment A8

**[flagged: plan amendment — A8, 2026-09-26; analysis run 2026-08-11] Fixed-width
EEG binning for the complexity–criticality models.** The pre-registered pipeline
(`step3b_branching_ratio_v2.py`) bins the population event series at the data-driven
mean inter-event interval, so dt ~ 1/event-rate. dt enters m through
m = exp(−dt/τ); event rate drives complexity (LZ vs events_per_epoch r = +0.938; LZ vs
dt r = −0.801; MPR_C vs dt r = −0.723). The binning parameter was on both sides of
the model, and within-subject label permutation cannot detect it. `step19` re-estimates
m at fixed dt = 5…11 samples (560/560 cells pass `mr_ok`); `step20` refits the models
with the variance control and 10,000 within-subject permutations per dt. This
amendment was made **after** the pre-registered results were examined, so it is
reported as a revision, with both specifications shown. Full reasoning is in
`docs_binning_confound_memo.md`.

**[result] Criterion 2 FAILS under fixed binning.** MPR_C ~ dist + variance:
p = 0.40–0.50 at all 7 dt, with 66–71% attenuation. Criterion 1 passes at all 7 dt
(p = 0.013–0.024, p_perm = 0.007–0.011). Criterion 3 limitation holds at all 7 dt.
Under pre-registered binning, criterion 2 had passed (p = 0.0102). Per the user's
decision (2026-09-26): **report both, lead with the failure.** The fixed-dt result is
primary; the data-driven result is reported as inflated by the shared binning path.
Scored by `step12_falsification_passfail.py` (binning-controlled arm now required;
the script exits if `fixed_dt_test_results.csv` is missing).

**[flagged: previously unreported result] step4 model 2.** Adding events_per_epoch
to the pre-registered model: LZ β_dist = +0.0752, p = 0.155; MPR_C β_dist = −0.0373,
p = 0.244. It was computed but never reported. It is now in manuscript §3.1 as the
statistical-control counterpart to the fixed-dt design control.

**[routine] Multi-dt verdict rule.** A criterion PASSES only if it passes at every
fixed dt, FAILS if it fails at every dt, and is MIXED otherwise. This is stricter than
the majority rule printed by `step20`. It does not affect any verdict (7/7 or 0/7 on
every criterion).

**[routine] Criterion 3 comparison metric (fixed-dt arm).** "LZ equally strong or
stronger" is judged on the scale-free within-subject |r|, with the 20% comparability
band from `step13`. **Flag:** the pre-registered arm's rule in `step12`
(`|β_lz| ≥ |β_mpr| or p_lz ≤ p_mpr`) compares raw β, which the manuscript (§2.5)
itself calls a scale artifact (LZ SD ≈ 3.5× MPR_C). It is left unchanged because it is
the registered rule, and its verdict (LIMITATION) matches the standardized comparison.

**[result] Manuscript redline.** `03_manuscript_A8_redline.docx` contains tracked
changes (author "Claude") to the abstract (still 250 words), Introduction, §2.3,
§2.5, §3.1, §3.2, Table 1, §3.3, Discussion, Limitations and the Figure 1 caption,
plus 3 comments flagging open issues: the title; the Maschke "replication" claim, which
conflicts with Friedman p = 0.68 for m across conditions; and the §3.6 m values, which
still use data-driven binning. `03_manuscript.docx` is unchanged.

**[owed] Follow-ups created by A8.** Resolve the Maschke claim. Retitle.
(step8 re-run and Figure 1: done 2026-09-27, see below.)
`PAPER_manuscript.md` is an older draft and has not been updated.

**[result, 2026-09-27] §3.6 m side re-run with fixed-dt m.** `step8_behavioral_split.py`
now takes the threshold as `argv[2]` and adds a fixed-dt m test (dt = 5…11; 10,000
group and 10,000 concentration-adjusted permutations per dt). It runs on its own RNG
stream (seed 20260927), so the original numbers reproduce exactly (checked at both
thresholds). Primary ≤20 split: m is lower in unresponsive subjects at 7/7 dt
(p_dir = 0.022–0.029), but survives concentration at **0/7** (p = 0.21–0.25; it was
0.045 under data-driven binning). Sensitivity ≤21 split: unadjusted p_dir = 0.056–0.068
and adjusted p = 0.52–0.57 at every dt. Complexity is independent of dt and unchanged.
The dissociation (complexity beyond dose, m not) now holds under **both** splits.
Outputs: `step8_fixed_dt_le20.csv`, `step8_fixed_dt_le21.csv`.
**Fix:** step8 printed a per-subject CSV path but never wrote the file; it now writes
it. The existing le20/le21 CSVs match the new output.

**[flagged: wording correction] "branching ratio tracks dose/concentration"** (Intro,
§3.6, Discussion) was never tested directly; it was inferred from m failing the dose
adjustment. It is reworded to "does not separate responsiveness once concentration is
controlled".

**[routine] Figure 1 regenerated from fixed-dt m.** Shown at dt = 8 samples (32 ms),
the middle of the 5–11 sweep (`FIG1_DT` in `make_figures.py`), with data-driven m as
open grey reference points. The caption notes that the pooled OLS line understates the
within-subject effect. The old figure is superseded.

**[routine] Complexity regression tests added.** `tests/test_complexity.py` pins LZ
and MPR_C to the CLAUDE.md reference values. The generating signals had not been
recorded; they were reconstructed and reproduce all six numbers to 3 decimals:
n = 5000, sine = 100 cycles over linspace, noise = N(0,1) from default_rng(0), chaos =
logistic map r = 4, x0 = 0.4. The tests also confirm that inline LZ76 equals antropy
(to 1e-12 relative), so results don't depend on whether antropy is installed. pytest
was added to requirements.txt.

---

## Culture arm scoping (2026-10-02)

**[flagged: dataset choice, user decision] Culture arm datasets.** The user chose
CRCNS hc-8 (Timme et al. 2016; rat hippocampal dissociated cultures, 435 recordings,
DIV 6–35, 60-electrode MEA, spike-sorted) as primary, and Zenodo 31085 (Charlesworth
et al. 2015, CC0, developing primary cultures with pharmacological perturbations) as
the second preparation, conditional on verifying bursting and regime spread. The
literature check found no prior work locating the MPR_C peak relative to m = 1 with
an estimator able to read m > 1. Must-cite precursors: Timme 2016 (Tononi complexity
peaks near criticality in these cultures; their model peak was at p ≈ 0.265 vs a
critical p ≈ 0.25); Lotfi et al. 2021 PRE 103:012415 (statistical complexity maximal
near a sync/desync transition in anesthetized rat cortex); Martin-Burgos et al. 2024
(hc-8 timescales only). No culture data has been downloaded or inspected.

**[flagged: hard rule 5 / manuscript correction] The transient-onset estimator's
validation used oracle onsets.** step10/step11 detect burst onsets on the full
simulated population A and only average the subsampled trace. Real data has only the
observed trace. `step21_transient_observed_validation.py` (T = 100k steps, the 7 step11
worlds, f = 0.02/0.05/0.10) compares the protocols:
- oracle: crossing 0.937–1.065; sign accuracy ≥ 0.875.
- observed-only onsets: crossing 1.039–1.275, undefined in 18 conditions, sign
  accuracy as low as 0.
- split-half (detect on one half, measure on the other): crossing 0.885–1.283,
  undefined in 29 conditions. **It does not fix the bias.**
Diagnostic: the true A at observed onsets is median 7–13 (vs 3 for the oracle), well
below saturation. So the failure is mostly false onsets (fluctuations that die out),
not latency into saturation. **Consequence:** manuscript §3.5 / Table 4 / Figure 4
("crosses at 1.011", "subsampling-robust by construction") overstate what the
estimator can do on real recordings and must be corrected. The estimator is NOT
deployable on hc-8 in its current form.

**[routine] Avalanche module and ground-truth tests.** `src/criticality/avalanches.py`
(avalanche extraction, discrete truncated power-law MLE with KS-chosen xmin, Shew
kappa, DCC) with `tests/test_avalanches.py`: zipf exponent recovered within 0.06;
critical Galton-Watson gives tau ≈ 1.5, alpha ≈ 2, DCC < 0.25, kappa ≈ 1; subcritical
kappa < 1. All 5 pass.

**[result, 2026-10-02] kappa / DCC under observed-only conditions (step22).** Same
bursting model and 7 worlds; f = 1, 0.10, 0.05, 0.02; bin = 1, 2, 4 model steps.
- kappa is **monotone in true m in every condition**, so it is a valid ordinal
  criticality axis within a fixed f and bin.
- The kappa = 1 crossing is **not** at m = 1. It ranges 0.82–1.29 or is undefined
  (base world: f = 0.02 gives kappa < 1 even at m = 1.3; f = 1 gives kappa > 1 even
  at m = 0.8). Absolute side cannot be read without knowing f.
- DCC argmin ranges m = 0.80–1.30. Uninformative about m in this model.
**Consequence:** none of the three estimators tested (stationary MR, transient onset,
kappa/DCC) gives absolute criticality from a subsampled bursting network with unknown
f. The offset relative to m = 1 is not identifiable in hc-8 with these methods. A
within-dataset ordinal design on kappa remains valid. Awaiting user decision on the
culture-arm framing.

---

## External review of the culture-arm draft (2026-10-03)

The author pasted a review of the A9 draft. All five points were checked against the
code and data, and all hold:
1. kappa rank reversals: across f (m = 0.95 at f = 0.10 gives kappa 0.879 vs m = 1.10
   at f = 0.02 gives 0.842), across worlds (5.5% of pairs at fixed f and bin), and
   4/84 curves not strictly monotone. The step22 check tolerated −0.01 decreases; it
   is now strict.
2. LZ with a median split on sparse counts shows a rate sawtooth on iid Poisson data
   (0.83 / 1.02 / 0.67 / 0.84 / 1.01 / 0.95 at means 0.3–3.0). A balanced split after
   tie-breaking gives 1.013–1.019 and still detects structure (0.922).
3. The two-lines test was mis-specified (quadratic vertex instead of Simonsohn's Robin
   Hood breakpoint; resampling rules unspecified).
4. The failure claim was overstated: it concerns these implementations on one model
   family with homogeneous thinning, not identifiability in hc-8.
5. The model check was not implemented; requirements.txt wrongly said MNE was
   unnecessary; README omitted steps 19–22.

**Actions:**
- A9 rewritten as r2 (still a draft): kappa reframed as an avalanche-distribution
  axis; LZ balanced split; full Robin Hood specification with breakpoint re-selection
  in every permutation and bootstrap, signed statistic T = min(z_b, −z_c), failed fits
  = −∞; non-rejection read as "no supported interior maximum"; a rank-stability gate
  (§9.1) and an expanded model check (§9.2) are required before lock.
- requirements.txt: mne==1.12.1 added; antropy marked optional.
- README: steps 19–22 and tests documented.
- step22: strict monotonicity, re-run (80/84 strictly monotone; all 4 exceptions are at
  f = 0.02).
- Manuscript redline: 13 further tracked changes correcting the transient-estimator
  claims (Abstract, Introduction, §2.7, §3.5, Discussion, Limitations, Conclusion,
  Figure 4 caption). The abstract is now 248 words.

**[correction] "sign accuracy as low as 0"** (entry of 2026-10-02 and A9 r1) was
technically true but misleading. Zeros occur only at onset threshold 3, where 2–5 of
the 8 m values have any estimate. At the usable threshold 2, the observed-onset
protocol's sign accuracy over defined estimates is 0.38–1.00, with up to 5 of 8 m
values misclassified at f = 0.02. The manuscript uses the latter figure.

**[note] Lost scratch files.** The temp-directory cleanup removed the earlier
`redline.py` and XML backups. The redline helper now lives at `tools/redline_docx.py`.

---

## Second review and kappa rank diagnostic (2026-10-03)

The second pasted review was checked; all six points hold. Confirmed directly:
balanced-split LZ gives 1.015 on a periodic one-spike-per-100-bins train (raw
0.0007) and 1.018 on a silent signal; the near-critical kappa ρ in step22 is 0.487
vs 0.776 pooled. Not raised by the reviewer: **the r2 MPR_C jitter has the same flaw
as the LZ fix.**

**Actions (A9 → r3, still a draft):**
- Tie handling withdrawn and made an open item, with deterministic candidates
  (surrogate-normalized LZ on integer counts, tie-aware ordinal patterns, smoothed
  rate) and a fixed selection rule gated by a tie/sparsity battery.
- Numerical fit failures get lower/upper p bounds, with an "unresolved" verdict if
  they straddle .05; "no interior candidate" stays a legitimate −∞.
- Freedman–Lane permutation within culture blocks, with DIV spline and Z as nuisance
  (Winkler et al. 2014).
- §9.1 fully specified before running.
- §9.2 acceptance: 1,000 datasets, rate ≤ .05 and one-sided Clopper–Pearson upper
  bound ≤ .065. B = 999 permutations is a logged deviation from the 10,000 rule
  (algorithm validation, not a significance claim), with a 25-dataset check at
  10,000.
- Described as a mixed-model adaptation of Simonsohn.
- §1 table reports both onset thresholds with coverage.
- Manuscript: threshold-3 coverage added inside the existing tracked insertion.

**[result] step23 kappa rank diagnostic:** all four pre-set criteria pass (overall ρ
0.860; tertiles 0.89–0.94; near ρ 0.640; reversal 6.7%). **However** ρ(kappa,
observed rate) = 0.991; partial ρ(kappa, m | observed rate) = 0.345; | population rate
= 0.114; and population rate alone ranks m better (0.928). In this model kappa is
close to a firing-rate proxy. **Flagged for the author:** the ordinal-criticality
framing and the C2 rate control may be incompatible. Proposed stage-1 rule: ρ(kappa,
rate) ≥ 0.90 in hc-8 means kappa is treated as not separable from rate.

---

## Third review (2026-10-03): modified option 1

**[flagged: author decision] Culture arm continues with limited stage-1 diagnostics**,
kappa described as an avalanche-distribution axis. The author relayed the third
review's recommendation as the choice. All five points checked and hold:
1. step23 observes Binomial(A, p_rec); neurons have no activity histories, so it is a
   variable-sampling-probability validation, not neuron-resolved.
2. My "rate control wipes out the analysis" claim is retracted: Pearson(kappa, rate)
   = 0.842, VIF = 3.44; kappa adds to rate for predicting m (t = 21.9 over linear
   rate, t = 6.2 over log rate, ΔR² = 0.027).
3. The tie battery conflated LZ and MPR_C. Criteria are now separate, and tie-aware
   ordinal patterns need their own definition and validation.
4. Freedman–Lane nuisance terms were inconsistent (DIV only in the reduced model).
   Now one nuisance set N is used everywhere. Whole-culture sign-flipping is primary,
   because detrending does not establish within-culture exchangeability.
5. **[bug fix]** step23's near-transition reversal used |Δm| ≥ 0.05 instead of the
   specified 0.10. Corrected: 9.3% (was reported as 16.2%). The four primary criteria
   are unchanged. step23 gains `--summary-only` to recompute from the saved CSV.
The ρ ≥ 0.90 stop rule is withdrawn and replaced by predetermined stage-1
diagnostics (A9 §8.1).

**[correction, 2026-10-03] I wrongly conceded the third review's point 2.** The VIF of
3.44 assumes a linear rate term, but kappa–rate is non-linear: a rate spline explains
90.5% of kappa's residual after linear rate control; Pearson against log rate = 0.984
(VIF 31.4); after a 4-df rate spline 2.8% of kappa variance remains, and kappa adds
ΔR² = 0.017 for m. Fix: in A9 §6.4, Z now enters as 3-df splines, consistent with
the §8.1 residual-kappa diagnostics. The real issue is a dilemma: linear control does
not remove rate, and flexible control leaves little kappa. The user asked (same day)
that feedback be checked and pushed back on where wrong, rather than accepted.

---

## Fourth review (2026-10-03, re-pasted 2026-10-22): spline convention and §8.1 definitions

1. **Partly upheld.** My "2.8% left / ΔR² 0.017" used patsy `cr(rate, df=4)` (3 model
   df), while A9 said "3 df" with no convention: a real mismatch. The reviewer's
   alternative (`cr(df=3)` = 2 model df: 11.1% left, ΔR² 0.105) under-fits the rate–m
   curve (R² 0.634 vs about 0.72), so its larger ΔR² mostly reflects uncaptured rate.
   Results plateau from 3 model df upward (left 0.028 / 0.025 / 0.024; ΔR² 0.017 /
   0.021 / 0.024 at 3 / 4 / 5 model df; log rate 0.032 / 0.027).
   **Fix:** A9 §6.1 now gives exact patsy bases: `cr(rate, df=5)` (4 model df),
   `cr(variance, df=5)`, `cr(DIV, df=4)`; "df" means the patsy argument and model
   df = df − 1. Sensitivity: df = 3, 4, 6 and log rate. Earlier log lines saying
   "4-df rate spline" or "3-df splines" are superseded by this convention.
2. **Upheld.** §8.1 rewritten with exact definitions:
   - residual kappa from one OLS fit with no culture term, then within-culture
     centering;
   - interleaved 60 s blocks for the split halves, which removes slow drift from the
     error term; adjacent-block dependence makes R an upper bound; 30 s and 120 s
     sensitivity;
   - reliability r_AB across recordings on centred residuals, Spearman–Brown R;
   - SE_meas = SD(ẽ_A − ẽ_B) / 2;
   - rule: C2 interpretable iff R ≥ 0.50 and SD_within ≥ 2 × SE_meas.
3. **Split lock** written into A9. Stage 1 freezes at a commit, pending author
   approval; stage 2 stays closed; stage-1-prompted changes are logged as
   data-informed amendments.

---

## 2026-11-24: manuscript corrections and tie-battery specification

**[flagged: not done] Stage-1 freeze commit.** The user said "go ahead and finish the
project". My attempt to record that as author approval of the A9 freeze was blocked
by the permission classifier, correctly: locking a pre-registration needs the author's
explicit approval. Nothing has been committed. A9 still reads "Pending: author
approval".

**[result] Section-2 manuscript fixes** (tracked changes, validated against the
original `03_manuscript.docx`; abstract 249 words). Every number was recomputed first:
- **Maschke:** "replication" claims removed in 3 places. The m Friedman test across
  conditions is p = 0.68 (data-driven) and 0.66–0.88 (each fixed dt 5–11). Maschke's
  contrast was anesthetic loss of consciousness; these data are sedation. Neither
  replication nor contradiction is claimed.
- **§3.1:** added the Friedman result and the LZ direction. LZ is higher than baseline
  in mild (+0.025, 17/20, p = 0.005), moderate (+0.027, p = 0.040) and recovery
  (+0.022, p = 0.003; Wilcoxon), the opposite of Schartner 2015. Conditions were
  recorded in a fixed order, so order and drug are confounded.
- **m < 1** is attributed to stationarity (MR slopes = lag autocorrelations), not to
  the exponential form, in 5 places (Abstract, §2.3, §3.5 ×3, Discussion).
- **Title:** not changed. That is the author's choice; the existing comment remains.

**[routine] Tie battery generators fixed in A9 §9.2.1 before any code.** See A9.

**[flagged: spec gap, logged before the battery ran] A9 §9.2.1 candidate (c)** is
listed for "both measures", but LZ on a continuous smoothed rate needs a binarization
rule the spec does not give. step24 therefore computes (c) for MPR_C only. If LZ
candidate (a) fails, the LZ contrast is dropped per A9 §4's selection rule, not rescued
by an improvised binarization.
**[routine] LZ phrase counts in step24** use antropy (normalize=False), verified
identical to inline lz76 on multi-symbol Poisson sequences (6/6 cases); inline is about
quadratic in length. The full battery (n = 360k, K = 20) runs about 1–3 h; the spec was
not reduced for compute.

---

## 2026-11-24 (later): Entropy submission package, plus corrections found while preparing it

**[flagged: user decision] Target journal Entropy (MDPI)**, methods framing; the culture
arm is not part of this paper. Guidelines were checked on the journal site: free-format
submission is accepted; abstract about 200 words; cover letter required; preprints
allowed. **APC is CHF 2600** on acceptance.

**[correction] LZ under sedation is not "opposite to the literature".** Boncompte et
al. 2021 (Br J Anaesth 127:254, same Chennu dataset) report an LZ *increase* at low
propofol doses in responsive subjects. Schartner 2015's decrease is from general
anaesthesia. The manuscript redline (§3.1) and the Entropy text are corrected; only the
persistence at recovery remains anomalous. This supersedes the memo's "LZ moves the
wrong way" framing.

**[correction] The Chennu data DOI in the old draft notes (10.17863/CAM.6875) is wrong**
(it resolves to a geology article). The correct DOI is **10.17863/CAM.68959** (CC BY 2.0
UK), confirmed on the repository item page and from the paper's data availability
statement.

**[correction] AR(1) censoring claim.** "Censored iff true τ > 48 bins" was not exact.
Re-run of step18_confirm: 4 of the 5 processes with true τ ≥ 50 bins were censored; one
(true 50, fitted 46.7) passed, because censoring acts on the fitted τ. All processes with
τ ≤ 40 passed. The window change altered m by ≤ 0.001; the time window recovers true m
within 0.007.

**[correction] Allen: "22 of 53 ≥ 0.98"** should be 22 of the 47 passing estimates
(3840 ms window, 4 ms bins).

**[note] Missing CSVs regenerated.** synthetic_mr_validation.csv (step9, seeded)
reproduces the fold table exactly. transient_onset_validation.csv (step10) is not
regenerated; the Entropy paper uses step21 instead.

**[note] The quality-filter m bound differs by arm** (EEG 0 < m < 1.2; spiking
0 < m < 1.5); the Entropy Methods state both.

**[note] allen_kmax_ceiling_diagnostic.csv** has a corrupted first header cell (a pasted
shell command prefixes "session_id"). It is harmless to the analyses, but should be
fixed in repo cleanup.

**[result] Submission files** (docs/entropy/):
- manuscript.md (source) and Vemula_Entropy_manuscript.docx, built by
  tools/build_manuscript_docx.py: abstract 200 words, 23 references (all cited, none
  missing, DOIs verified via PubMed/Crossref/DataCite), 4 figures, 3 tables; passes the
  OOXML validator;
- Vemula_Entropy_cover_letter.docx;
- figures from make_figures_entropy.py into figures_entropy/ (300 dpi PNG and PDF,
  palette validated for colour-blind safety).
Placeholders for the author: e-mail, funding, repository URL/DOI, AI-use wording, date.
