# H1

**H1: MISMATCH ON PHENOTYPE (RATE-BIMODAL ONLY) — model PD bimodal in 5/5 seeds (silhouette p<0.05), BIC k=2 in 5/5, but the split is PV+ (fast, regular) vs PV- (slower): model cells on the data's burst/pause centroids 0.00 ± 0.00 bursty vs data 20/45 = 0.44 (95% CI 0.31-0.59); own-cluster bursty fraction 0.73 ± 0.07; model 'bursty' cluster is PV- (Cramer's V 0.85 ± 0.13). Original rule (73 rows, no phenotype criterion): partial.**

Prediction written 2026-11-15T13:20:20 (git 63ec255d); compared 2026-11-18T13:41:41.

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

## (b) Bursty proportion and (c) phenotype

Model 0.73 ± 0.07 (own clusters, seed mean); model cells assigned to the nearest *data* cluster centroid: 0.00 ± 0.00. Data 20/45 = 0.44 (dedup 20/45). Model PV+ fraction 0.21.

## Model clusters vs PV type

Pooled cell-runs (rows: model cluster, columns: type):

| cluster | PV+ | PV- |
|---|---|---|
| bursty | 1 | 212 |
| regular | 59 | 18 |

Per seed: share of PV+ cells in the bursty cluster 0.00 ± 0.01, share of PV+ cells that are bursty 0.02 ± 0.04, Cramér's V 0.85 ± 0.13.

## (d) Per-cluster features vs data

Rank-sum tests, rat cluster (45 distinct units) vs the model cluster of the same name (model cells pooled over seeds):

| cluster | feature | data_median | model_median | n_data | n_model | mw_p |
|---|---|---|---|---|---|---|
| bursty | rate_hz | 21.8 | 15.2 | 20 | 213 | 0.00713 |
| bursty | cv | 2.01 | 0.866 | 20 | 213 | 1.49e-13 |
| bursty | lv | 0.437 | 0.896 | 20 | 213 | 4.21e-09 |
| bursty | frac_isi_lt10ms | 0.162 | 0.0166 | 20 | 213 | 0.000185 |
| bursty | short_isi_rate_relative | 0.587 | 0.42 | 20 | 213 | 8.96e-13 |
| bursty | pause_time_frac | 0.45 | 0.0317 | 20 | 213 | 1.44e-13 |
| bursty | burst_share | 0.267 | 0 | 20 | 213 | 6.04e-16 |
| bursty | log_isi_bimodality | 0.297 | 0.215 | 20 | 213 | 8.01e-09 |
| regular | rate_hz | 32.5 | 31.3 | 25 | 77 | 0.663 |
| regular | cv | 0.969 | 0.728 | 25 | 77 | 3.61e-05 |
| regular | lv | 0.342 | 0.501 | 25 | 77 | 0.049 |
| regular | frac_isi_lt10ms | 0.0895 | 0.0216 | 25 | 77 | 0.0202 |
| regular | short_isi_rate_relative | 0.203 | 0.305 | 25 | 77 | 0.00289 |
| regular | pause_time_frac | 0.122 | 0.00783 | 25 | 77 | 1.3e-11 |
| regular | burst_share | 0.00521 | 0 | 25 | 77 | 5.77e-09 |
| regular | log_isi_bimodality | 0.215 | 0.232 | 25 | 77 | 0.161 |

## Secondary: anesthesia-like slow oscillation

PD + slow oscillation (labelled secondary): seeds with p < 0.05: 5/5, BIC prefers k=2: 4/5, silhouette 0.41 ± 0.03, bursty fraction 0.76 ± 0.07, PV+ share of the bursty cluster 0.01 ± 0.02, Cramer's V 0.90 ± 0.13.

## Caveats

Rat: urethane-anesthetized 6-OHDA, ~100 s trains; model: 20 s per seed, unanesthetized PD regime in the primary test. Model seeds share one 58-cell network (only input noise differs). 28 of the 73 rat rows are copies of 28 other units (same title exported twice; the spike times of one copy are a subset of the other), so the effective n is 45. burst_share has length-dependent false positives (Poisson ≈ 0.10 at 300 s, 20 Hz). log_isi_bimodality follows the reference code (Pearson kurtosis in Sarle's denominator; values are about half of the usual coefficient).

