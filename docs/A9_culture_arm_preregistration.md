# Amendment A9 — Culture arm pre-registration (DRAFT r4, not locked)

Drafted 2026-10-02; revised 2026-10-03 (r2, r3) after two external reviews. Status:
**draft**.
It locks only after (a) the rank-stability gate (§9.1) and the full-pipeline model
check (§9.2) have run and their results are appended here, and (b) the author
approves. After lock, changes are logged in `docs/decision_log.md` as amendments.

**Lock status (split lock).** Stage 1 is frozen at the commit that adds this line,
before any data access. Frozen: §3 inclusion rules; §4 population spike train and
bin widths; §5 kappa definition; §6.1 spline bases and nuisance set; §8 stage
protocol; §8.1 diagnostic formulas, thresholds and the "not interpretable" outcome.
**Stage 2 stays closed** (§4 tie handling, §6.2–6.4 tests, §7, §9.2) until tie
handling is finalized and the full-pipeline validation passes. Any change prompted
by stage-1 results is logged as a **data-informed amendment**.
*(Pending: author approval of the freezing commit.)*

**Blinding statement.** As of r2, no CRCNS hc-8 or Zenodo 31085 file has been
downloaded, opened, or summarized. Only the public dataset descriptions have been
read (recording counts, DIV range, electrode count, mean firing rate 1.9 Hz).

**r2 changes (review of 2026-10-03):** kappa reframed as an avalanche-distribution
axis, with its observed rank reversals disclosed (§1). LZ binarization changed
because the median split creates rate artifacts on sparse counts (§4). The two-lines
test is now specified as Simonsohn's Robin Hood procedure in full, with resampling
rules (§6). The estimator-failure claim is narrowed to these implementations on this
model family (§1). A rank-stability gate is added and the model check is expanded
(§9).

**r3 changes (second review, 2026-10-03):**
- The r2 LZ fix is **withdrawn**. Random tie-breaking turns tie-dominated signals into
  noise: a periodic one-spike-per-100-bins signal goes from LZ 0.0007 to 1.015, and a
  silent signal scores 1.018. The same flaw applies to the r2 MPR_C jitter.
  Tie-handling is now an open item with deterministic candidates (§4).
- Numerical fit failures get conservative bounds instead of −∞ (§6.2).
- The adjusted permutation test uses Freedman–Lane with exchangeability blocks
  (§6.2, §6.4).
- The rank diagnostic is fully specified before running and reports strata and the
  near-transition region (§9.1).
- §9.2 acceptance rules are tightened, and the reduced-permutation deviation is
  stated (§9.2).
- Both onset thresholds are reported with coverage (§1).
- The test is described as a mixed-model adaptation of Simonsohn's procedure.

**r4 changes (third review, 2026-10-03):**
- step23 is relabelled as a variable-sampling-probability validation, not a
  neuron-resolved one (§9.1).
- The claim that the rate control makes C2 uninterpretable is retracted (VIF 3.44;
  kappa adds information beyond rate) (appendix).
- The pooled ρ ≥ 0.90 stop rule is replaced by predetermined stage-1 diagnostics
  (§8.1).
- The tie battery has separate LZ and MPR_C criteria, and any tie-aware ordinal
  measure needs its own definition and validation (§9.2.1).
- Nuisance terms are identical in observed and resampled fits (§6.1).
- Whole-culture sign-flipping is primary, with within-culture permutation as a
  sensitivity analysis (§6.2).
- The near-transition reversal metric is corrected to the specified spacing.

---

## 1. Scope: what this arm can and cannot claim

Validation on a saturating, depressing branching model with **homogeneous binomial
thinning** found that none of the three estimator implementations tested recovers the
critical point m = 1 under observed-only conditions:

| estimator | result under observed-only conditions | evidence |
|---|---|---|
| stationary MR | supercritical side folds below 1 | step9, manuscript §3.5 |
| transient onset | observed-onset crossing 1.04–1.28 (split-half 0.89–1.28). Sign accuracy over defined estimates, with coverage (defined estimates out of 8 m ≠ 1 per world): threshold 2 — f = 0.02: 0.38–0.71 (7–8/8); f = 0.05: 0.50–0.88 (8/8); f = 0.10: 0.62–1.00 (8/8). Threshold 3 — f = 0.02: 0.00–0.20 (4–5/8); f = 0.05: 0.33–0.71 (6–8/8); f = 0.10: 0.50–0.75 (8/8). Proposed cause: false onsets (supported by true A at onset = 7–13, not shown to be the sole cause) | step21 |
| kappa (Shew) | kappa = 1 crossing at 0.82–1.29 or undefined, depending on f and bin | step22 |
| DCC | argmin anywhere in m = 0.80–1.30 | step22 |

This does **not** show that absolute criticality is unidentifiable in hc-8. The
simulator does not represent fixed, spatially selected neurons with heterogeneous
firing, and other estimators exist. It does mean this arm makes **no claim about an
offset relative to m = 1**.

**kappa is an avalanche-distribution axis, not a validated criticality ordering.**
In step22 kappa increases with m within a fixed observation fraction f and parameter
world, but ranks reverse across those:
- Across f at the same bin: m = 0.95 at f = 0.10 gives kappa = 0.879; m = 1.10 at
  f = 0.02 gives 0.842.
- Across parameter worlds at fixed f and bin: 5.5% of (m_low, m_high) pairs reverse
  (871/15,876).
- 4 of 84 curves are not strictly increasing (the step22 check tolerated decreases of
  0.01; it is now strict).

hc-8 recordings plausibly differ in effective f and in network parameters, so fixed
bins alone do not make kappa a criticality ranking. Unless the rank-stability gate
(§9.1) passes, results are stated as "along the avalanche-distribution axis kappa".
Even if it passes, a peak along kappa does not establish criticality or offset
invariance.

## 2. Hypotheses

- **H_C1 (primary).** In hc-8, MPR_C has a supported interior maximum along kappa.
- **H_C2 (contrast).** LZ, under the §4 preprocessing, has no supported interior
  maximum along kappa. This is conditional on §9.2 showing that the preprocessing
  does not itself create LZ peaks from rate changes.
- **H_C3 (exploratory).** The DIV distribution of recordings in the MPR_C flat
  region. Descriptive only.

## 3. Data and inclusion (fixed before access)

- **Primary:** CRCNS hc-8 (Timme et al. 2016; doi:10.6080/K0PC308P), 435 recordings.
  Unit = recording; grouping = culture ID (random intercept).
- **Inclusion:** ≥ 20 sorted neurons; ≥ 30 min; ≥ 500 avalanches at the primary bin.
  Exclusions are reported with counts and reasons.
- **Second preparation (conditional):** Zenodo 31085 (Charlesworth et al. 2015, CC0).
  Included only if stage 1 finds ≥ 500 avalanches in ≥ 50% of recordings and a
  kappa IQR ≥ 0.10.

## 4. Signals (fixed bin widths in physical time)

- **Population spike train:** all included sorted neurons pooled.
- **Criticality bin:** 4 ms primary; 2 and 8 ms sensitivity; identical for every
  recording.
- **Complexity signal:** population count in 10 ms bins.
- **Ties (MPR_C and LZ): OPEN ITEM, must be resolved before §9.2.** Sparse counts
  are dominated by ties.
  - **Rejected:** random tie-breaking (r2) injects randomness. A periodic
    one-spike-per-100-bins signal and a silent signal both score LZ ≈ 1.02.
  - **Rejected:** the EEG-arm median split. On iid Poisson counts it gives a rate
    sawtooth: LZ = 0.83, 1.02, 0.67, 0.84, 1.01, 0.95 at means 0.3, 0.65, 0.75, 1.0,
    1.5, 3.0.
  - **Deterministic candidates to test:**
    - (a) LZ on the integer count sequence itself (no binarization), normalized by
      the mean LZ of within-recording shuffle surrogates, so that marginal
      rate/sparsity dependence cancels;
    - (b) MPR_C with ordinal patterns that encode ties as equalities (modified
      permutation entropy), or MPR_C normalized against shuffle surrogates;
    - (c) a Gaussian-smoothed rate (σ = 10 ms), which has no ties but imposes kernel
      structure.
  - **Selection rule (fixed now):** the first candidate, in order (a)/(b), then (c),
    that passes all tie/sparsity tests in §9.2.1 is adopted. If none passes, the LZ
    contrast (H_C2) is dropped and MPR_C is reported with the limitation stated.
    Silent or near-silent recordings (fewer than 500 avalanches) are already excluded
    (§3).
- **MPR_C:** d = 6, delay 1; same implementation as the EEG arm (pinned by
  `tests/test_complexity.py`), with the tie handling above.
- **LZ:** per the tie-handling item above. The EEG arm is unaffected, since continuous
  signals have no ties.

## 5. Axis

- **Primary:** kappa (`src/criticality/avalanches.kappa`, tau_ref = 1.5, 10 points)
  at 4 ms.
- **Descriptive only:** stationary MR τ (ms, A7 window rule); DCC.

## 6. Interior-maximum test: mixed-model adaptation of Simonsohn's (2018) two-lines test

Simonsohn's procedure (Robin Hood breakpoint) is specified for independent
observations with OLS. Here it is **adapted**: the interrupted regression has a
culture random intercept and permutation-based inference. It is therefore not an
exact implementation, and its error control is established by §9.2, not inherited
from the original paper. Applied identically to MPR_C and LZ. With covariates Z
(§6.4), every step uses them.

### 6.1 Procedure (one dataset)

**Spline bases (exact; patsy, in a model with an intercept):**
- rate: `cr(rate, df=5)` (4 model df);
- variance: `cr(variance, df=5)`;
- DIV: `cr(DIV, df=4)` (3 model df).

"df" below always means the patsy `df` argument. Model df = df − 1, because one
column is absorbed by the intercept. The rate flexibility is set where the step23
results stop changing (see appendix). Sensitivity: rate with `df` = 3, 4 and 6, and
`np.log(rate)`.

**Nuisance set N (identical everywhere):** `cr(DIV, df=4)`, plus Z when §6.4
applies. The same N enters the observed statistic, every resampled refit
and every bootstrap refit.

1. **Smoother.** Fit a cubic smoothing spline of y on kappa (pooled; smoothing chosen
   by GCV), with y residualized on N. ŷ_flat = fitted values within 1 SE of
   max ŷ. If max ŷ lies outside the 5th–95th kappa percentile, there is no interior
   candidate and the test **does not reject**.
2. **Interim breakpoint** x_c0 = median kappa within ŷ_flat. Fit the interrupted mixed
   model y = a + b·x_low + c·x_high + d·high + N·B + (1 | culture), where
   x_low = x − x_c if x < x_c (else 0), x_high = x − x_c if x ≥ x_c (else 0), and
   high = 1[x ≥ x_c]. Record z1 = |z_b| and z2 = |z_c|.
3. **Robin Hood breakpoint** x_c = the z2 / (z1 + z2) percentile of kappa within
   ŷ_flat. Refit the step-2 model at x_c.
4. **Decision.** Supported interior maximum if b > 0 and c < 0, each two-sided
   p < .05.

### 6.2 Test statistic and resampling

- **Signed statistic** T = min(z_b, −z_c). This is large and positive only when the
  lower line rises and the upper line falls; wrong-sign slopes make T negative.
- **Resampling scheme (hard rule 3): Freedman–Lane with culture blocks** (Winkler et
  al. 2014, NeuroImage 92:381). Shuffling kappa alone would break its relationship
  with the nuisance variables, so the null would not hold them fixed. Removing the
  age trend does **not** by itself make repeated recordings of one culture
  exchangeable, since residuals may be serially dependent across DIV.
  - **Primary: whole-culture sign-flipping.** Each draw multiplies the residual
    vector of each culture by an independent ±1 (y* = ŷ_N + s_c·e_c). This keeps the
    within-culture dependence structure intact. It assumes errors that are
    symmetric at the culture level, not exchangeable within culture.
  - **Sensitivity: within-culture permutation** of e (steps 1–3 below), which
    additionally assumes within-culture exchangeability after detrending.
  - Both schemes are validated in §9.2, including a scenario with serially
    correlated within-culture residuals. If they disagree on the data, the verdict
    follows the primary scheme, and the disagreement is reported.
  1. Fit the reduced model y ~ N + (1 | culture) by REML. N = nuisance: DIV as a
     natural cubic spline with 3 df, plus Z when §6.4 applies. Keep fitted values
     ŷ_N and residuals e.
  2. For each of 10,000 draws, sign-flip e by culture (primary), or permute e within
     culture (sensitivity).
  3. Form y* = ŷ_N + e_perm and rerun steps 1–3 **on y* with kappa unchanged**,
     including breakpoint selection. p_perm = (1 + #{T* ≥ T_obs}) / 10,001.
  - Cultures with a single recording form singleton blocks. They are sign-flipped
    under the primary scheme and left unpermuted under the sensitivity scheme.
- **Fit outcomes in any draw (observed or permuted):**
  - **No interior candidate** (step 1): a legitimate outcome meaning no U evidence.
    T = −∞.
  - **Numerical failure** (non-convergence, singular fit): T is unknown. p_perm is
    computed twice, with failures set to −∞ (lower bound) and to +∞ (upper,
    conservative bound). The verdict is reported only if both bounds fall on the same
    side of .05; otherwise it is **unresolved**. If the observed statistic itself is a
    numerical failure, the verdict is unresolved. Failure counts are always reported.
- **Bootstrap:** 10,000 culture-level resamples, rerunning steps 1–3. Reports the
  kappa range of ŷ_flat and the rejection rate. **The Robin Hood breakpoint is chosen
  to maximize power and is not interpreted as the peak location.** Location is
  reported only as the ŷ_flat range.

### 6.3 Interpretation

A rejection supports an interior maximum along kappa. A non-rejection means "no
supported interior maximum". It does **not** establish that the peak is absent or
unbracketed.

### 6.4 Confound control (hard rule 2)

Z = {population count variance, mean population firing rate}, separately and
together. **Each enters as a spline with the bases above, not linearly.** In step23,
kappa depends on rate non-linearly: after a linear rate term, 90.5% of the remaining
kappa variance is still explained by a rate spline. A linear control would leave the
rate confound inside kappa while appearing well-conditioned. Firing rate is the culture analogue of the A8 event-rate confound. Z
enters as nuisance in the Freedman–Lane reduced model (§6.2) and as covariates in
steps 1–3. The adjusted procedure is validated under firing-rate confounding in §9.2
before use.

## 7. Falsification criteria (culture arm)

1. **C1:** §6 rejects for MPR_C at 4 ms (parametric and permutation p < .05).
   Otherwise: no supported interior maximum.
2. **C2:** C1 still holds with Z (each covariate, and both). FAIL if not.
3. **C3:** if LZ also rejects, H_C2 fails and it is reported as a limitation. This is
   only interpretable if §9.2 passes for LZ.
4. **C4 (if Charlesworth is included):** rejection in both preparations. Location
   comparison is exploratory only.
5. **C5:** ≥ 15% of recordings pass inclusion.
6. **Robustness:** C1 is also evaluated at 2 and 8 ms. The claim rests on 4 ms;
   disagreement is reported, not resolved post hoc.

## 8. Two-stage blinding

1. **Stage 1:** inclusion statistics, kappa, MR τ, DCC, mean rate and variance only.
   Complexity code not run. Freeze (commit + hash). Decide on Charlesworth.
2. **Stage 2:** MPR_C and LZ, then §6–7 as written.

### 8.1 Predetermined stage-1 diagnostics (no complexity; replace the r3 ρ ≥ 0.90 rule)

Computed in stage 1, alongside kappa (bases as in §6.1):
1. **kappa–rate dependence:** Pearson and Spearman ρ(kappa, rate); R² of OLS
   `kappa ~ cr(rate, df=5)`; the same at df = 3, 4, 6 and for log rate.
2. **Residual kappa (definition).** e_r = kappa_r minus the fitted value from OLS
   `kappa ~ cr(rate, df=5) + cr(DIV, df=4)`, fit once on all included recordings
   with no culture term, so e retains culture differences.
   - Within-culture residual: ẽ_r = e_r − mean of e over the recordings of r's
     culture. Single-recording cultures are dropped from within-culture quantities
     and counted.
   - Reported: 1 − ICC of kappa (random-intercept model), SD(e), and SD_within =
     the pooled within-culture SD of ẽ.
3. **Neuron-count effect:** coefficient and culture-clustered SE of n_neurons in OLS
   `kappa ~ cr(rate, df=5) + cr(DIV, df=4) + n_neurons`.
4. **Interleaved split-half precision.**
   - Each recording is cut into consecutive 60 s blocks. Half A = odd blocks,
     half B = even blocks. Avalanches are extracted within blocks; those crossing a
     block boundary are dropped.
   - kappa_A, kappa_B, rate_A and rate_B are computed per half.
   - Half residuals e_A and e_B use the step-2 coefficients (from full recordings),
     evaluated at the half's own rate and the recording's DIV. Each is then centred
     within culture, giving ẽ_A and ẽ_B.
   - Reliability: r_AB = Pearson(ẽ_A, ẽ_B) across recordings; Spearman–Brown
     R = 2·r_AB / (1 + r_AB). The uncentred version is also reported.
   - Measurement SE of a full-recording residual: SE_meas = SD(ẽ_A − ẽ_B) / 2, which
     assumes independent half errors.
   - **Drift:** interleaving makes slow drift (slower than 60 s) common to both
     halves, so it is not counted as error. Faster non-stationarity, and dependence
     between adjacent blocks, are not removed. Adjacent-block dependence inflates
     r_AB, so R is treated as an upper bound. Sensitivity: 30 s and 120 s blocks.

**Rules (fixed now):**
- C2 is **interpretable** only if R ≥ 0.50 (60 s blocks) **and**
  SD_within ≥ 2 × SE_meas.
- Otherwise C2 is reported as "not interpretable: insufficient rate-independent kappa
  variation". That is an explicit outcome, neither a pass nor a failure.
- C1 is still run, and stated with the kappa–rate dependence alongside.
- No pooled-correlation cutoff is used to stop or to validate a criticality framing.
  Passing these diagnostics licenses only the avalanche-distribution description.

## 9. Validation before lock (no culture data)

### 9.1 Rank diagnostic for kappa (specification fixed 2026-10-03, before running)

A diagnostic, not a calibration. Passing supports kappa's ranking of m **under these
model assumptions only**. kappa is described as an avalanche-distribution axis
whatever the outcome.

- **Script:** `step23_kappa_rank_diagnostic.py`. Seed 20261003.
- **Recordings:** n_rec = 400 synthetic recordings, each with independent draws.
- **Dynamics:** step11 model (saturating, synaptically depressing branching process),
  population cap N = 2000. Per recording: m ~ U(0.80, 1.30), h ~ U(0.10, 0.40),
  tau_rec ~ U(150, 600) steps, u ~ U(0.03, 0.08) (the step11 ranges).
- **Time mapping:** 1 model step = 1 ms. Duration T = 300,000 steps (5 min) after a
  5,000-step burn-in. This is shorter than hc-8's ~60 min, so kappa is noisier here
  than it would be in the data (conservative).
- **Observation model:**
  - Every unit i in 1..N gets a fixed propensity w_i ~ LogNormal(0, 1), a skewed
    firing distribution.
  - Each recording observes a fixed random subset of n_obs ~ U{20, …, 142} units
    (the hc-8 range).
  - Each step's spikes are assigned to units in proportion to w, independently, so
    the observed count is Binomial(A(t), p_rec), with p_rec = Σ_obs w / Σ_all w.
  - Under mean-field dynamics this is exact. It produces heterogeneous effective
    observation fractions across recordings, which is the main source of the
    step22 reversals.
  - **Not modelled:** spatial structure, per-neuron refractoriness, heterogeneous
    burst participation. A pass says nothing about those.
  - **Scope (r4):** individual neurons have no activity histories; their weights only
    set p_rec. This is a variable-sampling-probability validation, **not** the
    neuron-resolved validation first proposed, which has not been done.
- **kappa:** `src.criticality.avalanches.kappa` at bins of 2, 4 and 8 ms; 4 ms is
  primary. Recordings with < 500 avalanches at a bin are excluded from that bin, and
  the exclusion count is reported.
- **Reported, at each bin:**
  1. overall Spearman ρ(kappa, m);
  2. ρ within tertiles of p_rec (sampling strata);
  3. ρ near the transition, |m − 1| ≤ 0.10;
  4. the reversal rate among pairs with |Δm| ≥ 0.10 (fraction where the higher-m
     recording has the lower kappa), overall and near the transition;
  5. partial Spearman ρ(kappa, m | p_rec, mean observed rate).
- **Diagnostic criteria (4 ms; all must hold):** overall ρ ≥ 0.80; every p_rec
  tertile ρ ≥ 0.70; near-transition ρ ≥ 0.50; overall reversal rate ≤ 10%.
- **Outcome use:**
  - Pass: the manuscript may say that kappa ordered m in a heterogeneous simulated
    population.
  - Fail: no ranking claim at all; peaks are described only along "the
    avalanche-distribution axis kappa", and the failure is reported.

### 9.2 Full-pipeline model check

Run §4–6 end to end (tie handling as resolved, spline, Robin Hood, Freedman–Lane
resampling) on
synthetic datasets with hc-8-like structure: about 60 cultures, repeated recordings,
heterogeneous f and rates.

| scenario | ground truth | requirement |
|---|---|---|
| monotone MPR_C along kappa | no interior max | acceptance rule below |
| known interior peak (pre-set effect size) | interior max | power reported; target ≥ 0.80 |
| firing rate drives both kappa and complexity, no true peak | no interior max | adjusted test (Z) rejection ≤ 0.05; unadjusted rate reported |
| sparse ties / rate changes only (no temporal structure) | no interior max in LZ or MPR_C | rejection ≤ 0.05 |
| sampling differences (f varies across recordings) | as above | as above |

**Status: specification incomplete. It must be completed and appended here before
§9.2 runs, and §9.2 runs only after the tie-handling item (§4) is resolved.** Still
to fix in writing before running:
- the signal generators for each scenario, including how a known MPR_C peak and a
  firing-rate confound are induced;
- the known-peak effect size, in SD units of MPR_C;
- the smoothing-uncertainty handling (GCV, with the 1-SE flat region propagated
  through resampling);
- an age-trended scenario that tests the within-culture exchangeability assumption
  (§6.2).

**9.2.1 Tie/sparsity battery (gates the §4 choice).** Each candidate is run on:
- silent and near-silent series;
- iid Poisson at means 0.1–5;
- a periodic sparse train (one spike per 50 / 100 / 200 bins);
- a rate-modulated Poisson;
- bursting-model output at f = 0.02–0.10.

**Generators and candidates (fixed 2026-11-24, before `step24_tie_battery.py` was written or run):**
- **Length:** n = 360,000 bins (1 h at the 10 ms complexity bin). Seed 20261124.
- **Silent:** all zeros. **Near-silent:** Poisson mean 0.001.
- **iid Poisson means:** 0.005, 0.01, 0.02, 0.1, 0.3, 0.65, 0.75, 1, 1.5, 3, 5. The low
  means match the periodic trains.
- **Periodic sparse:** one spike every P bins, P ∈ {50, 100, 200} (means 0.02, 0.01,
  0.005). Compared with iid at the same mean.
- **Rate-modulated Poisson:** λ(t) = μ(1 + 0.9 sin(2πt/200)), μ ∈ {0.3, 1, 3}.
  Compared with iid at the same μ.
- **Bursting model:** step11 base world, m = 1.0, 1 ms steps for 3.6 M steps,
  binomially observed at f ∈ {0.02, 0.05, 0.10}, summed into 10 ms bins. Compared
  with iid at the matched mean.
- **Shuffle surrogates:** K = 20 within-series random permutations.
- **Undefined:** fewer than 2 distinct values or fewer than 100 non-zero bins. Such
  series are flagged, never scored.
- **Candidate (a), LZ:** LZ_norm = LZ76(integer count sequence, no binarization) /
  mean LZ76 over the K shuffles. Null value 1.
- **Candidate (b), MPR_C:** ΔC = C(x) − mean C over the K shuffles. C is the existing
  MPR_C (d = 6, delay 1) with deterministic stable-sort ranking of ties. Null value 0.
  "Higher C" means ΔC > 0.02, i.e. beyond the iid tolerance band.
- **Candidate (c), both measures:** the same statistics on a Gaussian-smoothed rate
  (σ = 10 ms, on the 1 ms count series where one exists, otherwise on the 10 ms
  series). Tested only if (a) or (b) fails.
- Tie-aware ordinal patterns are **not** tested in this round, because they need
  their own definition first.

Criteria are **separate per measure**, because LZ and MPR_C respond to structure in
different directions (reference: sine MPR_C 0.211 > noise 0.028, whereas LZ is
sine < noise).
- **LZ candidates:**
  - iid within ±0.02 of the surrogate-null value at every rate;
  - periodic and rate-modulated series ≥ 0.2 below iid on the normalized scale;
  - silent/near-silent series flagged as undefined, never scored;
  - no periodic/iid reversal at any sparsity.
- **MPR_C candidates:**
  - iid within ±0.02 of its own surrogate-null value at every rate (no rate
    dependence on structureless data);
  - rate-modulated and bursting-model series differ from iid in the direction the
    continuous-signal references imply (higher C);
  - silent series flagged as undefined;
  - candidates are not required to rank periodic below iid.
- **Tie-aware ordinal patterns** change the pattern space (ordered set partitions
  instead of d! permutations) and therefore the normalization constants. Such a
  measure must be written out in full (pattern set, entropy normalization, Q0 for
  that pattern count) and get its own synthetic reference tests before the battery.
  It does not inherit `tests/test_complexity.py`.

**Acceptance rules (replacing "5% + 2 MC SE"):**
- 1,000 synthetic datasets per null scenario.
- A null scenario passes if the observed rejection rate is ≤ 0.05 **and** the
  one-sided 95% Clopper–Pearson upper bound is ≤ 0.065.
- This certifies the false-positive rate is below 6.5% at 95% confidence, not exactly
  5%. That is stated as such.
- Power for the known-peak scenario is reported with its 95% CI, and the target is
  ≥ 0.80.

**Permutation count (deviation from the standing 10,000-shuffle rule, logged):**
- Inside §9.2 each synthetic dataset uses B = 999 Freedman–Lane permutations, for
  compute.
- This validates the algorithm's error rate to within the Monte Carlo resolution of
  a 999-draw p-value. It makes no significance claim about data, which is what the
  10,000 rule governs.
- A random subset of 25 datasets per scenario is rerun with 10,000 permutations, and
  the B = 999 and B = 10,000 verdicts must agree in ≥ 24 of 25.
- The culture analysis itself always uses 10,000.

## 10. Out of scope

- Offset relative to m = 1.
- Any re-tuning after stage 2 begins.
- Mixing culture and in-vivo data in one model.
- Novelty framing: Timme et al. 2016 studied (Tononi) complexity versus criticality
  in these same cultures; Lotfi et al. 2021 showed statistical complexity is maximal
  near a sync/desync transition in rat cortex. The contribution here is limited to
  MPR_C with an LZ contrast, rate/variance controls and a valid interior-maximum
  test. It must be framed that way.

---

## Appendix — §9.1 result (run 2026-10-03, `step23_kappa_rank_diagnostic.py`, seed 20261003)

**All four pre-set criteria PASS at 4 ms** (400 recordings, none excluded):

| metric | value | criterion |
|---|---|---|
| overall ρ(kappa, m) | 0.860 | ≥ 0.80 |
| p_rec tertile ρ (low / mid / high) | 0.936 / 0.891 / 0.938 | ≥ 0.70 |
| near-transition ρ (\|m − 1\| ≤ 0.1, n = 161) | 0.640 | ≥ 0.50 |
| overall reversal rate | 6.7% | ≤ 10% |

Also: near-transition reversal 9.3% at |Δm| ≥ 0.10 as specified (an earlier printout
gave 16.2% because the code used 0.05; fixed). 2 ms and 8 ms are similar (overall ρ
0.856 / 0.847).

**But kappa is almost a monotone function of firing rate in this model**, which
undercuts the ordinal design:
- ρ(kappa, observed rate) = 0.991.
- Partial ρ(kappa, m | p_rec) = 0.923, but partial ρ(kappa, m | observed rate) = 0.345,
  and partial ρ(kappa, m | rate/p_rec) = 0.114.
- Population rate (observed rate / p_rec) ranks m **better** than kappa (ρ = 0.928 vs
  0.860).

**Correction (r4).** Two earlier statements — that kappa orders m "mainly because
both track firing rate", and that the C2 rate control would leave kappa "almost
without variance" — were overstated.
- **Correction to the r4 correction (2026-10-03).** r4 accepted the reviewer's VIF
  argument (Pearson 0.842, VIF 3.44 against linear rate). That argument does not
  hold:
  - The VIF is low only because rate enters linearly while kappa–rate is strongly
    non-linear (Pearson against log rate 0.984, VIF 31.4).
  - After a linear rate term, a rate spline explains 90.5% of the remaining kappa
    variance.
  - Results depend on spline flexibility until about 3 model df, then plateau:
    | rate basis | kappa variance left | ΔR² for m from kappa |
    |---|---|---|
    | linear | 0.291 | 0.316 |
    | cr(df=3) = 2 model df | 0.111 | 0.105 |
    | cr(df=4) = 3 model df | 0.028 | 0.017 |
    | cr(df=5) = 4 model df (chosen) | 0.025 | 0.021 |
    | cr(df=6) = 5 model df | 0.024 | 0.024 |
    | log rate | 0.032 | 0.027 |
    The cr(df=3) row under-fits the rate–m curve (R² 0.634 vs 0.72), so its larger
    ΔR² mostly reflects uncaptured rate, as with the linear row. An earlier r4 note
    called the 3-model-df fit "4-df" and the spec "3 df" without defining the
    convention; both are now defined above. The "t = 21.9 over linear rate" figure
    reflects the same non-linearity.
  - **The real issue is a dilemma:** linear control does not remove rate, and
    flexible control leaves little kappa variation. Under the step23 model, C2 is
    therefore expected to be weakly informative. The §8.1 reliability rule on
    spline-residual kappa decides this on the data.
- What the simulation shows is a strong, near-monotone kappa–rate association under
  these assumptions. It does not show that kappa reduces to rate in general. Rate
  can itself be a consequence of the network dynamics.
Whether this collinearity holds in cultures (where rate also reflects density and
excitability) can be checked in stage 1 without touching complexity.

**Status (r4):** the author chose to proceed with limited stage-1 diagnostics (§8.1).
The r3 candidate ρ ≥ 0.90 stop rule is withdrawn.
