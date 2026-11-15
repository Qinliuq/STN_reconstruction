# Held-out validation (H1–H4)

Experimental targets that are **not** used for fitting configuration D. The model is never re-tuned to these tests.

## Protocol
1. `priors.json`: fixes the input-regime parameters before any held-out run: seeds, durations, the slow-oscillation envelope, the D2/D1 routes and their magnitudes.
2. `model/run_all.py`: runs the simulations (specs in `model/specs/`) through `src/sweep.py`. It uses the same cfg writing and the same seed convention: `ctx_seed = s`, `gpe_seed = 100 + s`, with a fixed network layout. Output goes to `sweeps/out/heldout/`.
3. `model/predict_H1_H2.py`, `predict_H3.py`, `predict_H4.py`: compute the model metrics and write `results/heldout/H<n>/prediction.json`. Each file records the git hash and dirty flag, priors, the run configs, the spike-file hashes and the per-cell CSV hashes. Model code never imports `data_io.py`; `tests/test_isolation.py` checks this.
4. `compare/compare_H<n>.py`: refuses to run unless `prediction.json` exists and the files it lists still match their hashes. Only then does it load the data (`data_io.py`; set `HELDOUT_DATA_DIR` to the unpacked bundle). It applies a verdict rule written into each script's docstring before the comparison, and writes CSVs, `H<n>_summary.png`, `summary.md` and `comparison.json`.
5. `compare/make_readme.py`: builds `results/heldout/README.md`.

```bash
python -m pytest validation/heldout/tests -q
```

```bash
python -m validation.heldout.data_prep.verify_rat_clusters
```

```bash
python -m validation.heldout.model.run_all --jobs 4
```

```bash
python -m validation.heldout.model.predict_H1_H2 && python -m validation.heldout.model.predict_H3 && python -m validation.heldout.model.predict_H4
```

```bash
for t in 1 2 3 4; do python -m validation.heldout.compare.compare_H$t; done; python -m validation.heldout.compare.make_readme
```

## Shared metrics (`metrics/`)
One function serves both inputs (rat units and model cells).
- `spike_features.py`: exact port of `rat_6ohda/extract_features.py`, with the 8 features and the inclusion criteria (≥ 100 spikes, ≥ 99 positive ISIs). Recomputing from the extracted rat trains reproduces the saved features to 1e-13 (`results/heldout/H1/rat_verification.json`).
  - `log_isi_bimodality` keeps the reference definition, which puts the Pearson kurtosis in Sarle's denominator where Sarle uses the excess kurtosis. Its values are about half the usual coefficient: normal 1/6, uniform 0.21. The usual 0.555 threshold does not apply.
  - `burst_share` has false positives on long Poisson trains: about 0.10 at 20 Hz over 300 s. The rat "regular" cluster (0.0015) is therefore sub-Poisson.
- `clustering.py`: port of `cluster_test.py` (StandardScaler within the subset, Ward k=2, silhouette, 1000-draw Gaussian null with rng seed 0). The BIC check is a GaussianMixture on PCA-2 with full covariance and 10 inits, k = 1 vs 2.
  - Rat labels are reproduced exactly (ARI 1.0) and the silhouette is 0.367.
  - Null p95 is 0.321 and p is 0.005 here, against 0.316 and 0.002 in the reference. The null draws depend on the NumPy version.
  - BIC is k1 582.7 vs k2 555.7 (reference 583 vs 556).
- `phase.py`: zero-phase Butterworth band-pass and Hilbert phase, PPC0 (Vinck et al. 2010), preferred phase, Rayleigh p, and the slow-wave-fraction window rule.
- `response.py`: per-cycle window rates, per-cell Wilcoxon (p < 0.01), responder labels and trial-averaged profiles.
- `stats.py`: Mann–Whitney, KS, ECDF, Cramér's V / Fisher.

## Model input hooks (off by default)
`src/input_modulation.py` provides `cfg.ctx_rate_mod` and `cfg.gpe_rate_mod`; both default to `None`. With the hooks off, a default run is bit-identical to the code before the hooks: spike times and GPe rates were compared on a PD seed-1 run.
- `slow_osc`: a ~1 Hz up/down envelope on the cortical arrival rate, normalized to mean 1.
- `blocks`: a pre/on/post protocol. For the cortex it multiplies the arrival rate; for the GPe it multiplies the baseline (target) rate of each group, so the STN→GPe feedback stays active.

The saved spikes now also include `ctx_rate` (cortical population arrival rate, 1 ms), `gpe_beta_phase` / `gpe_in_burst` (per GPe step) and, when the hooks are on, `ctx_env` / `ctx_so_phase`.

## Data notes
- **Rat duplicates.** In `unit_cluster_SSonly.csv`, 28 (file, label) titles appear twice. One copy's spike times are always a subset of the other's (same unit exported twice), so the 73 rows correspond to 45 distinct units. `extract_rat_mat.py` writes both copies under one label. `data_io.rat_ss_trains()` merges coincident spikes (< 1 µs) to recover the larger copy.
  - The 73-row analysis is reported as published. A 45-unit sensitivity analysis gives a 20/25 split, silhouette 0.341, p = 0.047, and BIC still prefers k = 2.
- **No activated-cortex epochs.** The rat unit trains are 100 s epochs that lie entirely in slow-wave-dominated EEG. Activated windows (slow fraction < 0.4) exist only outside them, so beta locking in activated cortex cannot be tested. Beta over all epochs is reported as a secondary result.
- **No EEG for two sessions.** pjx286c and pjx289c have no EEG channels; their 8 units are excluded from H3.
- **ChR2 D1 count.** `binned_response.csv` holds 50 D1 units (27 excited), not 60.
- **ChR2 time course.** The KS4 D2 trial-aligned profiles are in the bundle only as figures, so time courses are compared with the binned early / late / post windows.
