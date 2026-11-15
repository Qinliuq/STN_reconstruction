# H2

**H2 (pooled check): MISMATCH — similar distributions for 0/3 features; medians model vs rat: rate_hz 17.7 vs 25.7; cv 0.837 vs 1.32; burst_share 0 vs 0.0327.**

Prediction written 2026-11-15T13:20:20 (git 63ec255d); compared 2026-11-15T13:20:29.

## Tests

| feature | data | model_cells | n_a | n_b | median_a | median_b | mw_U | mw_p | ks_D | ks_p |
|---|---|---|---|---|---|---|---|---|---|---|
| rate_hz | rat_73 | all | 290 | 73 | 17.7 | 25.7 | 6.62e+03 | 7.62e-07 | 0.36 | 2.97e-07 |
| rate_hz | rat_73 | PV+ | 60 | 73 | 32.4 | 25.7 | 2.82e+03 | 0.00448 | 0.545 | 1.47e-09 |
| rate_hz | rat_73 | PV- | 230 | 73 | 15.7 | 25.7 | 3.8e+03 | 1.92e-12 | 0.523 | 1.71e-14 |
| rate_hz | rat_dedup_45 | all | 290 | 45 | 17.7 | 29 | 3.68e+03 | 2.45e-06 | 0.387 | 9.56e-06 |
| rate_hz | rat_dedup_45 | PV+ | 60 | 45 | 32.4 | 29 | 1.57e+03 | 0.159 | 0.472 | 1.1e-05 |
| rate_hz | rat_dedup_45 | PV- | 230 | 45 | 15.7 | 29 | 2.11e+03 | 3.28e-10 | 0.563 | 1.13e-11 |
| cv | rat_73 | all | 290 | 73 | 0.837 | 1.32 | 3.66e+03 | 5.45e-18 | 0.702 | 1.34e-28 |
| cv | rat_73 | PV+ | 60 | 73 | 0.733 | 1.32 | 547 | 1.11e-13 | 0.802 | 8.64e-22 |
| cv | rat_73 | PV- | 230 | 73 | 0.861 | 1.32 | 3.11e+03 | 5.44e-16 | 0.696 | 1.55e-26 |
| cv | rat_dedup_45 | all | 290 | 45 | 0.837 | 1.32 | 2.39e+03 | 8.31e-12 | 0.687 | 1.86e-18 |
| cv | rat_dedup_45 | PV+ | 60 | 45 | 0.733 | 1.32 | 348 | 8.88e-11 | 0.789 | 2.2e-16 |
| cv | rat_dedup_45 | PV- | 230 | 45 | 0.861 | 1.32 | 2.05e+03 | 1.44e-10 | 0.681 | 1.93e-17 |
| burst_share | rat_73 | all | 290 | 73 | 0 | 0.0327 | 4.44e+03 | 1.06e-18 | 0.491 | 2.1e-13 |
| burst_share | rat_73 | PV+ | 60 | 73 | 0 | 0.0327 | 723 | 1.24e-12 | 0.6 | 1.33e-11 |
| burst_share | rat_73 | PV- | 230 | 73 | 0 | 0.0327 | 3.71e+03 | 1.03e-15 | 0.465 | 2.1e-11 |
| burst_share | rat_dedup_45 | all | 290 | 45 | 0 | 0.0327 | 2.31e+03 | 3.23e-16 | 0.569 | 2.64e-12 |
| burst_share | rat_dedup_45 | PV+ | 60 | 45 | 0 | 0.0327 | 346 | 1.35e-12 | 0.678 | 9.72e-12 |
| burst_share | rat_dedup_45 | PV- | 230 | 45 | 0 | 0.0327 | 1.97e+03 | 8.49e-14 | 0.54 | 1.06e-10 |

## Note

Pooled check only: the model was not tuned to these distributions and will not be. Model cells are 5 seeds x 58 cells of one network (not independent); rat units are urethane-anesthetized.

