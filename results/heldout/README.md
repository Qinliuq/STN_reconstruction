# Held-out validation of configuration D

Experimental targets that were never used for calibration. Each model prediction (`H<n>/prediction.json`, with git hash, configuration and spike-file hashes) was written before the comparison script read the data; the comparison scripts refuse to run without it. Input-regime parameters were fixed in `validation/heldout/priors.json` before any run. Configuration-D parameters were not changed. Code: `validation/heldout/` (see its README).

| test | data source | model metric | model | data | statistic | verdict | caveats |
|---|---|---|---|---|---|---|---|
| H1 | Rat 6-OHDA single units (unit_cluster_SSonly.csv, 73 SS rows / 45 distinct units) | bimodality (silhouette vs Gaussian null, BIC), bursty fraction, cluster vs PV+ | bimodal 5/5 seeds (BIC k2 5/5); bursty 0.73 ± 0.07; V(PV+) 0.85 ± 0.13 | bimodal (sil 0.367, p 0.005; BIC k2); bursty 33/73 = 0.45 | silhouette vs Gaussian null, GMM BIC, Wilson CI, Cramer V | **partial** | urethane anesthesia; 28/73 rat rows are copies of other units (n_eff 45); model seeds share one 58-cell network; 20 s model trains vs ~100 s rat trains |
| H2 | Rat 6-OHDA single units (73 SS rows) | pooled rate / CV / burst share (pooled check) | rate_hz 17.7; cv 0.837; burst_share 0 | rate_hz 25.7; cv 1.32; burst_share 0.0327 | Mann-Whitney + KS (medians shown) | **mismatch** | pooled check only, never tuned; anesthesia; species/state of rates differ |
| H3 | Rat 6-OHDA units + IpsiEEG (37 units with EEG) | PPC + preferred phase, slow 0.5-2 Hz and beta 13-30 Hz | slow PPC 0.0862; beta PPC 0.00337 (ctx) / 0.00367 (GPe) | slow PPC 0.0806; beta (all epochs) -4.62e-05; activated: 0 units | Mann-Whitney on PPC (medians shown) | **match** | model slow oscillation is an imposed input envelope (priors.json); EEG polarity vs cortical firing unknown; no rat activated epochs inside unit epochs (beta test on activated cortex not possible) |
| H4 | Mouse D1/D2-MSN ChR2 (binned_response.csv; D2 wave_clus) | responder sign / fraction / magnitude, time course, PV enrichment | pending | | | pending | mouse, dopamine-intact, anesthetized, multi-unit (~1-4 Hz); compressed 1/2/2 s protocol (slow recovery / after-effects untested); route magnitudes a priori |

One-line verdicts:

- H1: PARTIAL — model PD bimodal in 5/5 seeds (silhouette p<0.05), BIC k=2 in 5/5; bursty fraction 0.73 ± 0.07 vs data 33/73 = 0.45 (95% CI 0.34-0.57); bursty cluster vs PV+: Cramer's V 0.85 ± 0.13.
- H2 (pooled check): MISMATCH — similar distributions for 0/3 features; medians model vs rat: rate_hz 17.7 vs 25.7; cv 0.837 vs 1.32; burst_share 0 vs 0.0327.
- H3: MATCH — slow PPC match (model median 0.0862 vs rat 0.0806); beta PPC to cortical input not testable (model 0.00337 vs rat activated-epoch nan, 0 rat units qualify); beta to GPe drive not testable (model 0.00367). Secondary, slow-wave epochs: beta mismatch (model PD+SO 0.0043 vs rat all-epoch -4.62e-05).

Details: `H<n>/summary.md`, figures `H<n>/H<n>_summary.png`, tables `H<n>/*.csv`.
