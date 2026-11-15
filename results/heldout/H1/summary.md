# H1

**H1: PARTIAL — model PD bimodal in 5/5 seeds (silhouette p<0.05), BIC k=2 in 5/5; bursty fraction 0.73 ± 0.07 vs data 33/73 = 0.45 (95% CI 0.34-0.57); bursty cluster vs PV+: Cramer's V 0.85 ± 0.13.**

Prediction written 2026-11-15T13:20:20 (git 63ec255d); compared 2026-11-15T13:20:27.

## (a) Bimodality

Per seed (58 cells; seeds share the network):

| seed | n | silhouette | null_p95 | p | bic_k1 | bic_k2 | bursty_frac | frac_PVp_in_bursty | cramers_V | fisher_p |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 58 | 0.353 | 0.349 | 0.048 | 469 | 434 | 0.621 | 0 | 0.653 | 7.25e-07 |
| 2 | 58 | 0.407 | 0.368 | 0.009 | 478 | 435 | 0.759 | 0 | 0.905 | 1.02e-10 |
| 3 | 58 | 0.408 | 0.374 | 0.013 | 468 | 450 | 0.707 | 0 | 0.793 | 6.94e-09 |
| 4 | 58 | 0.417 | 0.365 | 0.002 | 480 | 427 | 0.793 | 0.0217 | 0.895 | 6.2e-10 |
| 5 | 58 | 0.424 | 0.352 | 0.003 | 469 | 431 | 0.793 | 0 | 1 | 1.12e-12 |

Mean ± SD over seeds: silhouette 0.40 ± 0.03, null p95 0.36 ± 0.01, p 0.015 ± 0.019, BIC k1 472.6 ± 5.8 / k2 435.5 ± 8.6. Pooled over seeds (290 cell-runs, not independent): silhouette 0.422, p 0.000, BIC k1 2285.6 / k2 1982.4.

Data (port of cluster_test.py, labels reproduced exactly): silhouette 0.367, null p95 0.321, p 0.005 here (reference 0.316 / 0.002; the null draws depend on the NumPy version); BIC k1 582.7 / k2 555.7 (reference 583 / 556). With the 28 duplicated unit copies removed (45 units): silhouette 0.341, p 0.047, BIC k1 365.4 / k2 356.5, 20/25 split.

## (b) Bursty proportion

Model 0.73 ± 0.07 (own clusters, seed mean); model cells assigned to the nearest *data* cluster centroid: 0.00 ± 0.00. Data 33/73 = 0.45 (dedup 20/45). Model PV+ fraction 0.21.

## (c) Bursty cluster vs PV+

Pooled cell-runs (rows: model cluster, columns: type):

| cluster | PV+ | PV- |
|---|---|---|
| bursty | 1 | 212 |
| regular | 59 | 18 |

Per seed: share of PV+ cells in the bursty cluster 0.00 ± 0.01, share of PV+ cells that are bursty 0.02 ± 0.04, Cramér's V 0.85 ± 0.13.

## (d) Per-cluster features vs data

Rank-sum tests, rat cluster vs the model cluster of the same name (model cells pooled over seeds):

| cluster | feature | data_median | model_median | n_data | n_model | mw_p |
|---|---|---|---|---|---|---|
| bursty | rate_hz | 21.4 | 15.2 | 33 | 213 | 0.0078 |
| bursty | cv | 2 | 0.866 | 33 | 213 | 2.5e-20 |
| bursty | lv | 0.405 | 0.896 | 33 | 213 | 2.26e-13 |
| bursty | frac_isi_lt10ms | 0.0461 | 0.0166 | 33 | 213 | 0.00013 |
| bursty | short_isi_rate_relative | 0.591 | 0.42 | 33 | 213 | 1.71e-19 |
| bursty | pause_time_frac | 0.451 | 0.0317 | 33 | 213 | 2.4e-20 |
| bursty | burst_share | 0.294 | 0 | 33 | 213 | 7.5e-23 |
| bursty | log_isi_bimodality | 0.324 | 0.215 | 33 | 213 | 1.78e-15 |
| regular | rate_hz | 31.1 | 31.3 | 40 | 77 | 0.631 |
| regular | cv | 0.995 | 0.728 | 40 | 77 | 8.34e-07 |
| regular | lv | 0.28 | 0.501 | 40 | 77 | 6.04e-06 |
| regular | frac_isi_lt10ms | 0.0188 | 0.0216 | 40 | 77 | 0.947 |
| regular | short_isi_rate_relative | 0.162 | 0.305 | 40 | 77 | 3.93e-05 |
| regular | pause_time_frac | 0.132 | 0.00783 | 40 | 77 | 5.77e-15 |
| regular | burst_share | 0.0015 | 0 | 40 | 77 | 1.64e-07 |
| regular | log_isi_bimodality | 0.248 | 0.232 | 40 | 77 | 0.276 |

## Secondary: anesthesia-like slow oscillation

PD + slow oscillation (labelled secondary): seeds with p < 0.05: 5/5, BIC prefers k=2: 4/5, silhouette 0.41 ± 0.03, bursty fraction 0.76 ± 0.07, PV+ share of the bursty cluster 0.01 ± 0.02, Cramer's V 0.90 ± 0.13.

## Caveats

Rat: urethane-anesthetized 6-OHDA, ~100 s trains; model: 20 s per seed, unanesthetized PD regime in the primary test. Model seeds share one 58-cell network (only input noise differs). 28 of the 73 rat rows are copies of 28 other units (same title exported twice; the spike times of one copy are a subset of the other), so the effective n is 45. burst_share has length-dependent false positives (Poisson ≈ 0.10 at 300 s, 20 Hz). log_isi_bimodality follows the reference code (Pearson kurtosis in Sarle's denominator; values are about half of the usual coefficient).

