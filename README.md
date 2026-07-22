# Complexity and criticality in cortical dynamics — analysis code

Analysis code for the manuscript *"Complexity tracks approach to criticality but the
offset is not measurable in accessible cortical data: a boundary result and a fit-window
artifact in branching-ratio estimation."*

The study asks whether the offset between the statistical-complexity peak and the critical
point (branching ratio m = 1) is invariant across species, states, and recording modalities.
The main results are a boundary result (accessible cortex is subcritical and does not bracket
the peak), a structural limit of the stationary MR estimator, and a methodological finding
that a bin-specified fit-quality filter can manufacture an apparent subcritical ceiling.

## Data (not included here)

Both datasets are openly available from their original providers and are **not** redistributed
in this repository. Download them and note the paths.

- **Human EEG (propofol sedation):** Chennu et al. (2016), University of Cambridge Repository.
  EEGLAB `.set` files plus `datainfo.mat`. Place them in one folder; the scripts take that
  folder as their first argument.
- **Mouse spiking (Neuropixels):** Allen Institute Visual Coding Neuropixels dataset. Retrieved
  by `allen.py` into a local cache; the branching-ratio scripts read the cached NWB/HDF5 files
  directly with h5py.

## Environment

Python 3 (developed on Homebrew-managed Python, macOS/Apple Silicon). Install dependencies:

    pip install -r requirements.txt

See `requirements.txt` for the package list and a note on pinning the exact versions you ran.

## How the scripts fit together

Scripts take the data folder as their first command-line argument and write their outputs (CSV)
into that same folder, so later scripts find earlier outputs. Run them from, or point them at,
one working directory.

### EEG arm

1. `step2_batch_extract.py <data_dir>` — computes Lempel-Ziv (LZ76, implemented inline) and
   Rosso-Plastino statistical complexity (MPR_C) per channel/epoch.
   Writes `complexity_epochs.csv` and `complexity_summary.csv`.
2. `step2b_variance.py <data_dir>` — raw per-epoch signal variance (the amplitude confound
   control). Writes `variance_summary.csv`.
3. `step3b_branching_ratio_v2.py <data_dir>` — branching ratio via the Wilting-Priesemann MR
   estimator. Writes `branching_summary_v2.csv`.
   (`step3_branching_ratio.py` is the superseded v1; see "Superseded scripts" below.)
4. `step4_join_and_test.py <data_dir>` — joins complexity to branching ratio, fits the
   linear mixed-effects models, produces the §3.1 rising-limb result.
   Reads `complexity_summary.csv` and `branching_summary_v2.csv`.
5. `step5_variance_confound.py <data_dir>` — re-runs the complexity-criticality models with
   variance included (the §3.2 confound control).
6. `step6_permutation.py <data_dir>` — within-subject label-permutation significance
   (10,000 shuffles).
7. `step13_standardized_betas.py <data_dir>` — z-scored cross-measure comparison (Table 1).
8. `step17_eeg_kmax_sensitivity.py <data_dir>` — EEG fit-window sensitivity sweep; confirms the
   EEG arm is not affected by the fit-filter artifact. Writes `eeg_kmax_sensitivity.csv`,
   `eeg_kmax_summary.csv`, `branching_summary_A7.csv`.

### Behavioural (responsiveness) split — §3.6, exploratory

9. `step8_behavioral_split.py <data_dir>` — primary responsiveness split (correct-count rule).
10. `step8b_chennu_classification.py <data_dir>` — reproduces the original Chennu drowsy/responsive
    classification as the sensitivity variant. Writes `step8b_chennu_classification.csv`.

### Mouse spiking arm

11. `allen.py` — downloads/caches the Allen Neuropixels sessions.
12. `step7_allen_mscan_h5py.py <cache_dir>` — branching-ratio scan of the spiking data via the
    direct h5py loader. (`step7_allen_mscan.py` is the earlier AllenSDK-based version; the h5py
    version is the one used, because the SDK's NWB reader is unusable on the current stack.)
13. `step14b_locomotion_scope_v2.py <cache_dir>` — running-speed quantification bounding the
    §3.4 quiescence claim. (`step14_locomotion_scope.py` is the superseded v1.)
14. `step15_allen_mscan_stationary.py <cache_dir>` — the time-specified-window branching scan.
    Writes `allen_m_scan_stationary.csv` (+ residual-check CSV).
15. `step16_kmax_ceiling_diagnostic.py <cache_dir>` — the fit-filter diagnostic: shows the
    bin-specified KMAX filter imposes the apparent ceiling. Writes
    `allen_kmax_ceiling_diagnostic.csv` and `allen_kmax_ceiling_summary.csv`.
    (This is the core of the §3.3 methodological result.)

### Synthetic validation — §2.7, §3.5, Figures 3–4

16. `step9_synthetic_validation.py` — MR estimator against known generative m; demonstrates the
    non-injective fold (Figure 3). (Standalone; no data input.)
    *(Note: this file is named `step9_synthetic_validation.py`; some documentation refers to it
    as `step9_synthetic_mr_validation.py`. Same script.)*
17. `step10_transient_onset_prototype.py` — the transient-onset estimator that recovers the sign
    of (m − 1) (Figure 4). Standalone.
18. `step11_param_robustness_sweep.py` — parameter-robustness sweep for the transient-onset
    estimator. Standalone.
19. `step18_filter_artifact_demo.py` and `step18_confirm.py` — AR(1) ground-truth demonstration
    that the censoring is a filter artifact while the magnitude bias is a data property. Both
    are self-contained (the generator is inlined).

### Falsification scoring and figures

20. `step12_falsification_passfail.py <data_dir>` — scores the pre-registered falsification
    criteria as explicit pass/fail.
21. `make_figures.py` — regenerates all four figures. Reads `branching_summary_v2.csv`,
    `complexity_summary.csv`, the synthetic-validation CSVs, and (for Figure 2 panel B)
    `allen_m_scan_stationary.csv` and `allen_kmax_ceiling_diagnostic.csv`. **Set the `PROJ` and
    `OUT` paths at the top to `.` (or your data folder) before running.** Figure 2 prints its
    data source; "computed from allen_*.csv" confirms it read the real data rather than the
    fallback values.

## Superseded scripts

Two scripts are earlier versions kept for provenance. The paper uses the v2 in each case.
Do not run the v1 versions for the reported numbers.

- `step3_branching_ratio.py` — superseded by `step3b_branching_ratio_v2.py`.
- `step14_locomotion_scope.py` — superseded by `step14b_locomotion_scope_v2.py`.
- `step7_allen_mscan.py` — the AllenSDK-based scan, superseded by `step7_allen_mscan_h5py.py`.

## Reproducing the figures only

If you just want the figures from the summary CSVs (without re-running the full pipeline),
you need these files in one folder: `branching_summary_v2.csv`, `complexity_summary.csv`,
`synthetic_mr_validation.csv` (or your `step9` output), `transient_onset_validation.csv`,
`allen_m_scan_stationary.csv`, `allen_kmax_ceiling_diagnostic.csv`. Then set `PROJ`/`OUT` in
`make_figures.py` and run it.

## Citation

If you use this code, please cite the manuscript (citation to be added on publication).
