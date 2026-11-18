# H4

**H4: PARTIAL — D1_ctx match, D2_ctx match, D2_gpe mismatch. Data: D2 binned 37/69 inhibited (0 excited), responders x0.56; D1 27/50 excited, x1.79.**

Prediction written 2026-11-18T11:01:14 (git 5cced658); compared 2026-11-18T11:16:51.

## Comparison

| scenario | route | pathway | verdict | model_pop_ratio | data_pop_ratio | model_frac_inh | model_frac_exc | data_inh_exc | model_resp_ratio | data_resp_ratio | model_onset_late_post | data_early_late_post | resp_frac_PVp | resp_frac_PVm | type_V | type_fisher_p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| D1_ctx | D1-MSN, textbook route: cortical/thalamic drive increased | D1 | match | 1.17 ± 0.01 | 1.43 | 0.00 ± 0.00 | 0.62 ± 0.07 | 0/50 inh, 27/50 exc | 1.21 ± 0.01 | 1.79 | 1.20 / 1.19 / 1.00 | 2.08 / 1.79 / 1.34 | 0.58 ± 0.13 | 0.63 ± 0.06 | 0.04 | 0.459 |
| D2_ctx | D2-MSN, alternative route: cortical/thalamic drive reduced | D2 | match | 0.87 ± 0.01 | 0.68 | 0.45 ± 0.15 | 0.00 ± 0.00 | 37/69 inh, 0/69 exc | 0.83 ± 0.02 | 0.56 | 0.86 / 0.86 / 1.01 | 0.51 / 0.62 / 0.92 | 0.50 ± 0.10 | 0.44 ± 0.17 | 0.05 | 0.467 |
| D2_gpe | D2-MSN, textbook route: GPe rate reduced | D2 | mismatch | 1.19 ± 0.01 | 0.68 | 0.00 ± 0.00 | 0.71 ± 0.03 | 37/69 inh, 0/69 exc | 1.23 ± 0.02 | 0.56 | 1.20 / 1.20 / 1.00 | 0.51 / 0.62 / 0.92 | 0.80 ± 0.14 | 0.68 ± 0.04 | 0.10 | 0.0815 |

## Data targets

- D1: source=binned, n=50, n_exc=27, n_inh=0, frac_resp=0.54, ci=(np.float64(0.40398634328388083), np.float64(0.6703056539821737)), median_ratio_resp=1.792, early=2.078, late=1.794, post=1.342, pop_ratio=1.426
- D2: source=binned, n=69, n_exc=0, n_inh=37, frac_resp=0.536, ci=(np.float64(0.4197820076036184), np.float64(0.6488600870236277)), median_ratio_resp=0.561, early=0.511, late=0.618, post=0.922, pop_ratio=0.679
- D2_wc: source=wave_clus (crosstalk removed), n=75, n_exc=3, n_inh=25, frac_resp=0.333, median_ratio_resp=0.474

## Notes

binned_response.csv holds 50 D1 units (the task text says 60). The KS4 trial-aligned profiles are available only as figures in the bundle, so the time course is compared with the binned windows (early / late / post; window definitions in the ChR2 analysis) and the model protocol is compressed (1 s pre / 2 s on / 2 s post): the 20 s recovery and 10-20 s after-effects are not tested. Route magnitudes are a priori (priors.json), not fitted.

## Caveats

Mouse, dopamine-intact, anesthetized, multi-unit clusters (low rates ~1-4 Hz) vs a rat-based healthy model state without anesthesia; striatal stimulation reaches the STN through routes the model represents only as input-rate changes.

