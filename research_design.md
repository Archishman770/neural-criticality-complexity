# Research Design — Algorithmic and Statistical Complexity Around Neural Criticality

**Status:** Pre-registered as of the date this file is locked in project knowledge.
Amendments must be logged in `docs/decision_log.md`, not made silently by editing
this file in place — if a change is needed, add an "Amendment" entry below rather
than rewriting the original text.

## 1. Research question

Is the offset between the **statistical-complexity peak** (MPR_C) and the **true
branching-ratio critical point** (m_MR = 1) an *invariant* — conserved across
species, brain states, and recording modalities?

The peak-at-criticality phenomenon itself is established (closest precursor:
PMC10664178). This project does not re-derive that peak. The contribution under
test is whether the *offset* — its sign, magnitude, and/or existence — is
conserved when the same pipeline is applied across independent axes of variation
(species, brain state, recording modality). PMC10664178 does not test
cross-axis invariance; that gap is what this design is built to close.

## 2. Primary hypothesis

**H1:** Complexity (MPR_C) tracks distance from criticality
(`dist = |1 - m_MR|`) in a consistent direction and with a consistent offset,
across datasets that vary independently in species, brain state, and recording
modality.

**H0 (null):** No such relationship, or a relationship whose sign/magnitude is
not conserved across the tested axes (i.e., it is dataset-specific rather than a
general property of neural dynamics near criticality).

LZ is tracked in parallel as a secondary, contrasting measure — since LZ is
monotonic in disorder rather than peaked, it is **not** expected to show the same
relationship, and a null result for LZ where MPR_C shows a real effect is itself
informative (it argues the effect is about structure, not raw signal
randomness/amplitude).

## 3. Measures

| Axis | Measure | Role |
|---|---|---|
| Complexity (structure) | MPR_C (Bandt-Pompe ordinal patterns, Martín-Platino-Rosso statistical complexity) | Primary complexity measure; non-monotonic, peaks between order and randomness |
| Complexity (randomness) | Normalized Lempel-Ziv (LZ) | Contrast measure; monotonic in disorder, not expected to peak |
| Criticality (ground truth) | Branching ratio via Wilting-Priesemann multistep-regression (MR) estimator | Subsampling-robust criticality estimate; m_MR = 1 is the critical point |
| Confound | Raw signal variance | Amplitude-only control, computed independently of both complexity and event-detection pipelines |

Avalanche/event detection for the branching ratio uses **Option A**
(standard-deviation thresholding on z-scored signal, |z| > threshold peaks). This
is a fixed default, not re-litigated per dataset; Option B (negative-peaks-only)
is an explicit ablation if ever run, never a silent substitution.

## 4. Datasets

**Primary (in hand):** Chennu propofol sedation dataset — 20 subjects, 4
conditions (baseline / mild / moderate / recovery sedation), 91-channel EEG,
250 Hz, pre-epoched. Provides the brain-*state* axis (sedation level) within a
single species and modality.

**Planned extensions** (each adds an independent axis to the invariance test,
using the same pipeline and the same hard rules below):
- Polysomnography sleep EEG — additional brain-state axis (sleep stages) within
  the same modality (EEG), different paradigm than sedation.
- Rodent LFP / spiking data — species axis, and (for spiking) a modality axis
  (spikes vs. field potentials).
- Psychedelic EEG (optional) — a third brain-state axis, same modality.

The invariance claim is only as strong as the number of genuinely independent
axes tested. A result that holds only within the Chennu dataset supports "complexity
tracks criticality in this dataset" — not the invariance hypothesis. Do not present
single-dataset results as if they answer the core research question.

## 5. Falsification criteria

Stated in advance, per Section 3 hard-rules discipline in `CLAUDE.md`. For **each
dataset**, the pipeline's result is checked explicitly against these and reported
as **pass/fail**, not just as an effect size:

1. **Primary relationship (per dataset):** In the mixed-effects model
   `mpr_c ~ dist + (1|subject)`, restricted to quality-passing criticality cells
   (`mr_ok == True`), the coefficient on `dist` must be in the hypothesized
   direction (negative — complexity higher as distance-from-criticality
   shrinks) at p < .05. **Fails** if the sign is wrong or the effect is not
   significant.
2. **Survives the confound control:** The same relationship must remain
   significant (p < .05 on `dist`) after adding raw variance as a covariate
   (`mpr_c ~ dist + variance + (1|subject)`). **Fails** if `dist` loses
   significance once variance is included — this would indicate the effect is
   an amplitude artifact, not criticality-specific.
3. **LZ contrast:** LZ is not expected to show the same relationship. This is
   not a fail condition for H1, but if LZ shows an equally strong or stronger
   relationship than MPR_C, that undermines the claim that structure
   (not randomness) is doing the work, and must be reported as a limitation, not
   omitted.
4. **Cross-dataset invariance (only assessable once ≥2 independent-axis datasets
   are in):** The sign of the `dist` effect must match across datasets. If the
   *offset itself* (not just the sign) is the claim under test, the magnitude
   must fall within a pre-specified tolerance band (to be fixed before the second
   dataset's results are unblinded — do not set the tolerance after seeing the
   data).
5. **Quality gate:** Any dataset where fewer than ~15% of subject×condition cells
   pass the branching-ratio quality filter (`mr_ok`) is reported as
   underpowered/inconclusive for that dataset, not silently excluded from the
   writeup.

If a dataset fails any of 1–3, that is a genuine falsification result for that
dataset and must be reported as such — not re-analyzed with a different model
specification until it passes (see Section 7, pre-registration discipline).

## 6. Variance confound control (mandatory, every time)

Never report a complexity–criticality correlation without this alongside it.
Three tests, all pre-specified (mirrors `step5_variance_confound.py`):

- **A.** `complexity ~ dist + variance + (1|subject)` — does `dist` survive once
  variance is in the model?
- **B.** `variance ~ dist + (1|subject)` — does variance itself track
  criticality? (If yes, and complexity doesn't beat variance in C, the complexity
  measure may be doing no independent work.)
- **C.** Head-to-head within-subject correlation: criticality vs. complexity vs.
  criticality vs. variance. Complexity should outperform variance as a predictor
  for the result to be interpreted as criticality-specific rather than an
  amplitude effect.

Read: if `dist` dies in (A), variance is significant in (A), and variance tracks
criticality as well as complexity does in (C) — the effect is an amplitude
confound, not evidence for H1. If `dist` survives (A), that's the real signal.

## 7. Statistical analysis plan (pre-registered)

Run on quality-passing cells only (`mr_ok == True`), joined on
(subject, condition):

1. **Primary:** `mpr_c ~ dist + (1|subject)` and `lz ~ dist + (1|subject)`
   (mixed-effects, random intercept per subject).
2. **Event-rate control:** `complexity ~ dist + events_per_epoch + (1|subject)`
   — does the relationship survive controlling for raw event rate (a detection
   artifact risk flagged during branching-ratio development)?
3. **Within-subject:** per-subject Spearman correlation of complexity vs. m_MR
   across conditions (only for subjects with >2 distinct m_MR values), then a
   Wilcoxon signed-rank test on the distribution of subject-level rhos against
   zero.
4. **Between-subject:** Spearman correlation of subject-mean complexity vs.
   subject-mean m_MR.
5. **Variance confound control (Section 6):** run alongside every instance of
   1–4 that reports a complexity–criticality correlation.

**Amendments log:** any change to this plan after it is locked in project
knowledge gets appended below with a date and one-line rationale — the original
plan is never edited in place.

> _Amendments:_
> - _A8 (2026-09-26): EEG complexity–criticality models re-run at fixed bin width (dt = 5–11 samples) because data-driven binning put event rate on both sides of the model; made after pre-registered results were seen; both specs reported, fixed-dt primary; criterion 2 fails under it. See `docs/decision_log.md`. (A1–A7 are cited elsewhere but not yet entered here.)_

## 8. Permutation invariance testing

For any cross-condition or cross-group comparison reported as significant:
**minimum 10,000 label shuffles.** Procedure: shuffle condition (or group) labels
within subject, recompute the test statistic (e.g., the `dist` coefficient or the
between-group difference), build the null distribution, and report the
permutation p-value alongside the parametric one. If runtime is a constraint,
reduce scope (fewer channels or epochs per shuffle) — never reduce the shuffle
count below 10,000. If this constraint is ever at risk of being cut for
expedience, that must be flagged explicitly in output, not silently skipped.

## 9. Branching-ratio quality filter (fixed, set before results are inspected)

A subject×condition cell is `mr_ok = True` only if:
- The MR exponential fit converges (no exception),
- Fitted tau ≥ `MIN_TAU_BINS` (2.0 bins) — rules out degenerate near-zero-lag
  fits,
- Fitted tau ≤ `MAX_TAU_FRAC` (0.8) × `KMAX_BINS` — rules out fits that didn't
  decay within the fit window,
- `0 < m_MR < 1.2` — rules out non-physical values.

Cells failing any criterion are marked `mr_ok = False` with a `reason` string and
dropped from all correlational analyses (Sections 5–7). They are **not** dropped
silently — the drop rate and reasons are reported per dataset (Section 5,
criterion 5).

## 10. Reporting discipline

For every dataset run through this pipeline, the writeup states, in order:
1. Pass/fail on each falsification criterion (Section 5).
2. The confound-control result (Section 6), not just the primary correlation.
3. The permutation p-value alongside the parametric one, for any comparison
   claimed as significant.
4. The `mr_ok` drop rate and reasons.

A result is only described as supporting the invariance hypothesis (H1) once
it has been checked against Section 5's criteria across at least two
independent-axis datasets (Section 4) with a pre-specified tolerance for the
offset. Single-dataset results are reported as exactly that — evidence toward,
not confirmation of, the core claim.
