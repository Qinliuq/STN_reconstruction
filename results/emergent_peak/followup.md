# Emergent beta peak: POST-HOC follow-ups FU1 / FU2, 2026-11-23

Requested after the Step 2 results (`summary.md`); everything here is **post hoc**. Same seeds (1–5, fixed network),
duration (3 s), criteria and analysis code (`analyze.py`, unchanged; `analyze_followup.py`). Spec
`sweeps/emergent_peak_followup.py` and the analysis script were committed before any follow-up result existed.
50 runs, all succeeded. Files: `followup_conditions.csv`, `followup_profiles.csv`, `followup_profiles.png`,
`pv_pair_ratio_all.csv`.

| condition | definition |
|---|---|
| FU1 | configuration D (graded drive centre 0.65 / width 0.15 / floor 0.2, σ 0.1, configuration-D receptors) with the measured *Pvalb* PV+ placement |
| FU2_F | F (measured placement, *Gria1*, *Gabra1* + per-type ratio, uniform drive) with configuration D's GPe→STN topography (σ 0.1, even homes) |
| FU2_M0/1/2 | M (graded drive at 0.572 / 0.634 / 0.697, drive width 0.219) with configuration D's GPe→STN topography; healthy reference = FU2_F healthy (healthy beta depth is uniform) |
| Dref, Dref_uni | configuration D unchanged (healthy; PD graded; PD uniform drive), reruns for the PV+/PV− table |

**Reproducibility:** Dref reproduces the committed configuration-D numbers (PD 0.65/0.25 ratio 2.37, uniform drive
1.02; rates within 0.0001% of `results/gaba_ratio`). FU1 healthy is the same configuration as Step 2's ref healthy and
gives identical output.

## Verdicts (same criteria; post hoc)

- **FU1 — Table 3's "PV+ gradient necessary", original form: the peak survives with the measured placement, weaker.**
  Mean peak 0.70, 0.65/0.25 ratio 1.65 (≥ 1.5), 3/5 seeds in 0.50–0.80 → criterion not met on the seed count.
  Configuration D itself (Dref: 0.65, 2.37, **also 3/5 seeds**) fails the seed count too, so the per-seed rule is
  strict for this network. Measured vs exaggerated placement: ratio 1.65 vs 2.37, peak position unchanged.
  - So under the graded drive the measured *Pvalb* gradient keeps the peak at 0.70 and above 1.5, but contributes
    less contrast than the configuration-D cluster. Whether it is *necessary* needs the graded-drive + uniform-PV+
    rerun in configuration D (the Table 3 row, 1.25, predates the per-type ratio); not run.
- **FU2_F — uniform drive, measured receptors, configuration-D topography: no peak.** Mean peak 0.10, ratio 1.02
  (= Dref_uni 1.02), 2/5 seeds. Restoring the topography does not let the measured distributions localize beta under a
  uniform drive; the Step 2 negative for F does not depend on σ 0.22.
- **FU2_M (assumes motor-territory origin of beta): the peak appears at the two lower loci.**

  | locus | mean peak | per-seed peaks | 0.65/0.25 | seeds in range | criterion |
  |---|---|---|---|---|---|
  | M0, 0.572 | 0.65 | 0.55 0.65 0.65 0.70 0.60 | 2.45 | 5 | **met** |
  | M1, 0.634 | 0.65 | 0.55 0.75 0.60 0.75 0.65 | 2.25 | 5 | **met** |
  | M2, 0.697 | 0.70 | 0.45 0.75 0.80 0.65 0.95 | 1.98 | 3 | not met (seed count) |

  With the data topography (Step 2 M) the same drives gave ratios 0.54–0.79. So the data GPe→STN topography (σ 0.22,
  homes 0.41–0.98) and/or the rate drop that came with it is what removed the peak. FU2 changes σ, the home mapping and
  the rates together; this run cannot separate them.
- **Rates:** all follow-up conditions are within 4% of configuration D (largest: FU2_F healthy −3.9%). This confirms
  that the 17–25% rate drop in Step 2 came from the data GPe topography, not from the receptor profiles or placement.

## Per-condition table (PD; healthy rows in `followup_conditions.csv`)

| condition | mean peak | per-seed peaks | 0.65/0.25 | seeds in [0.50, 0.80] | PD − healthy max | beta pow | sync | bursts | rates PV+ / PV− (Hz) | rate change vs D |
|---|---|---|---|---|---|---|---|---|---|---|
| FU1 | 0.70 | 0.70 0.05 0.60 0.40 0.70 | 1.65 | 3 | 0.70 | 0.031 | 0.017 | 0.16 | 33.8 / 16.1 | −0.4% |
| FU2_F | 0.10 | 0.65 0.65 0.05 0.85 0.05 | 1.02 | 2 | 0.65 | 0.069 | 0.033 | 0.18 | 34.8 / 15.9 | −0.1% |
| FU2_M0 | 0.65 | 0.55 0.65 0.65 0.70 0.60 | 2.45 | 5 | 0.65 | 0.045 | 0.021 | 0.18 | 34.5 / 15.8 | −1.0% |
| FU2_M1 | 0.65 | 0.55 0.75 0.60 0.75 0.65 | 2.25 | 5 | 0.65 | 0.048 | 0.025 | 0.17 | 34.3 / 15.5 | −2.2% |
| FU2_M2 | 0.70 | 0.45 0.75 0.80 0.65 0.95 | 1.98 | 3 | 0.70 | 0.039 | 0.021 | 0.17 | 34.2 / 15.6 | −1.8% |
| Dref | 0.65 | 0.25 0.65 0.95 0.65 0.60 | 2.37 | 3 | 0.70 | 0.035 | 0.016 | 0.16 | 32.9 / 16.4 | 0.0% |
| Dref_uni | 0.50 | 0.05 0.55 0.50 0.05 0.45 | 1.02 | 2 | 0.05 | 0.070 | 0.029 | 0.17 | 33.2 / 16.7 | +1.3% |

Observation (not a criterion): in FU2_F the absolute beta peaks ventrally but the PD − healthy contrast peaks at 0.65,
because its healthy profile (measured receptors) is lowest in the middle.

## PV+ / PV− pairwise shared beta, all conditions so far (`pv_pair_ratio_all.csv`)

Per seed: PV+-pair beta cross-power / PV−-pair beta cross-power. Mean ± SD over seeds as asked, plus robust readouts:
pairwise cross-power can be near zero or negative (healthy runs; one abs PD seed has PV− = −0.2 ×10⁻⁶), which makes
the per-seed ratio unstable there.

| set | condition | state | ratio, mean ± SD of seeds | ratio of seed means | PV+ − PV− (×10⁻⁶), mean ± SD | seeds PV+ > PV− |
|---|---|---|---|---|---|---|
| Step 2 | F | PD | 2.71 ± 0.58 | 2.74 | 15.7 ± 7.3 | 5 |
| Step 2 | F_cfgD | PD | 2.78 ± 1.13 | 2.86 | 12.6 ± 10.5 | 4 |
| Step 2 | noPV | PD | 2.01 ± 0.74 | 1.98 | 8.8 ± 7.4 | 4 |
| Step 2 | noRec | PD | 3.02 ± 0.85 | 2.84 | 14.9 ± 6.9 | 5 |
| Step 2 | abs | PD | −2.67 ± 13.9 (one PV− < 0) | 3.96 | 10.5 ± 7.7 | 5 |
| Step 2 | sig | PD | 2.83 ± 0.45 | 2.87 | 13.5 ± 5.6 | 5 |
| Step 2 | ref | PD | 2.13 ± 0.87 | 1.85 | 10.6 ± 5.3 | 5 |
| Step 2 | M0 | PD | 2.84 ± 1.66 | 2.80 | 9.4 ± 6.1 | 4 |
| Step 2 | M1 | PD | 3.21 ± 1.91 | 2.74 | 9.1 ± 5.5 | 5 |
| Step 2 | M2 | PD | 2.35 ± 0.48 | 2.39 | 9.5 ± 5.0 | 5 |
| follow-up | FU1 | PD | 4.36 ± 2.65 | 3.29 | 10.8 ± 4.0 | 5 |
| follow-up | FU2_F | PD | 1.99 ± 0.33 | 1.93 | 13.5 ± 5.2 | 5 |
| follow-up | FU2_M0 | PD | 4.04 ± 1.56 | 3.77 | 20.0 ± 7.6 | 5 |
| follow-up | FU2_M1 | PD | 3.08 ± 0.47 | 3.10 | 17.8 ± 6.7 | 5 |
| follow-up | FU2_M2 | PD | 2.25 ± 1.67 | 2.43 | 9.7 ± 12.2 | 4 |
| follow-up | Dref | PD | 3.27 ± 1.47 | 3.35 | 12.1 ± 8.6 | 5 |
| follow-up | Dref_uni | PD | 2.38 ± 0.34 | 2.40 | 20.1 ± 10.2 | 5 |
| Step 2 | F | healthy | 17.5 ± 35.6 | 2.88 | 3.2 ± 6.6 | 3 |
| Step 2 | F_cfgD | healthy | 1.88 ± 1.09 | 1.91 | 1.7 ± 2.1 | 4 |
| Step 2 | noPV | healthy | 1.92 ± 3.82 | 2.21 | 1.6 ± 4.3 | 3 |
| Step 2 | noRec | healthy | 9.47 ± 12.9 | 3.73 | 4.4 ± 3.0 | 5 |
| Step 2 | abs | healthy | 4.07 ± 3.01 | 2.87 | 2.7 ± 2.8 | 4 |
| Step 2 | sig | healthy | 1.41 ± 1.43 | 1.36 | 0.6 ± 2.2 | 3 |
| Step 2 | ref (= FU1 healthy) | healthy | 1.71 ± 1.80 | 1.66 | 1.6 ± 3.4 | 3 |
| follow-up | FU2_F | healthy | 1.37 ± 1.06 | 1.37 | 1.0 ± 2.9 | 3 |
| follow-up | Dref | healthy | 4.93 ± 4.69 | 2.61 | 4.1 ± 3.4 | 4 |

**PD (17 conditions, 85 seed-runs):** PV+ pairs carry more shared beta than PV− pairs in 81 / 85 seed-runs and in
every condition (≥ 4 / 5 seeds); ratio of seed means 1.85–3.96 (configuration D 3.35). This holds whatever the drive
topography (graded, uniform, motor locus), PV+ placement (cluster, measured, uniform), receptor profiles and GPe
topography. Even with uniformly placed PV+ (noPV) the ratio is 1.98. It is the most robust result of this test.
**Healthy (9 distinct conditions):** PV+ > PV− in 32 / 45 seed-runs (FU1 healthy = ref healthy counted once); per-seed ratios are unreliable because PV−
pair beta is near zero; do not quote healthy ratios as mean ± SD.

## What this changes in `summary.md` / the results-draft note

- Table 3 "PV+ gradient necessary": in its original (graded-drive) form the peak position survives the measured
  placement (0.70) with lower contrast (1.65 vs 2.37); necessity itself is still untested (needs graded drive +
  uniform PV+ in configuration D).
- Step 2's M negative is due to the data GPe→STN topography (σ 0.22 is an upper bound) and/or its rate drop: with
  configuration D's topography a drive at the digitized somatomotor locus (0.57–0.63) gives the peak with the
  measured receptors and placement (criterion met, 5/5 seeds). This is still a drive-imposed peak, not an emergent one.
- The uniform-drive negative (F) is robust to the topography.
- Main robust result: PV+ > PV− pairwise shared beta (ratio ≈ 2–4) in every PD condition tested.
