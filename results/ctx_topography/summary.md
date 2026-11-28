# Cortical→STN topography test, 2026-11-28

## Setup

**Switch.** `cfg.ctx_topo_variant` (default `None`, which leaves configuration D unchanged). Code is in `src/ctx_topography.py`.
- Every scale is normalised to mean 1 over the 58 cells.
- The same variant is used in healthy and PD.

**Variants.**

| variant | cortical arrival rate | cortical conductance |
|---|---|---|
| T0 | 0.65 − 0.5x (current rule) | uniform |
| T1 | uniform | uniform |
| T2 | uniform | 0.568 − 0.273x (re-digitised Jeon 2022 Fig. 7H.2 Gria1 fit) |
| T3 | max(T_motor(x), 0.3) | uniform |

- T_motor is the Plantinga motor territory from `stn_drive_topography`. The 0.3 floor is a literal floor, applied before normalisation.

**Fixes and crossing.**
- F1 (`loc_grid 'centred'`) and F2 (`gpe_coverage_norm`) are on throughout.
- Each variant is crossed with pallidal depth:
  - graded: configuration D default;
  - uniform: `beta_width 1e9`, as in `Dref_uni_pd`.

**Runs.**
- Input seeds 1–5 (ctx s, gpe 100 + s); one network (loc 1, conn 4); 3 s per run.
- One healthy reference per variant. Healthy GPe depth is 0.1 everywhere, so that one reference matches both PD depths.
- Spec: `sweeps/ctx_topography.py`. Analysis: `analyze.py`. Outputs: `rates.csv`, `summary.csv`, `profiles.json`, `maps.png`.

**Check.** With `ctx_topo_variant = 'T0'`, the healthy seed-1 spikes are bit-identical to `R3_H_gradient_healthy_s1`, which used the old code path.

**Measures** (as in `results/edge_artefacts/rate_normalised.md`):
- rate-normalised pairwise beta, C_ij/(ν_iν_j), kernel 0.12;
- across-episode PPC₁ magnitude, pooled over seeds;
- ν·m = kernel-mean rate × √(Δf × normalised beta), with Δf = 1.95 Hz. This is the shared beta-band rate modulation of the local population, in Hz.
- Seed-noise check: peak − profile mean, compared with the SEM across the 5 seeds at the peak.
- A bin is flagged when n_eff < 8. Here that is x = 0.05 and 0.95 (n_eff 7.9) and every per-type PV+ bin except 0.65.

## Operating point (stop rule: mean PV+ or PV− rate within 10% of T0)

| | PV+ Hz | PV− Hz | PV+ vs T0 | PV− vs T0 |
|---|---|---|---|---|
| T0 healthy / graded PD / uniform PD | 37.7 / 32.7 / 33.4 | 20.3 / 16.8 / 17.0 | — | — |
| T1 healthy / graded / uniform | 40.3 / 35.4 / 35.5 | 19.5 / 16.1 / 16.2 | +7.1 / +8.4 / +6.4 % | −3.7 / −4.2 / −4.7 % |
| T2 healthy / graded / uniform | 39.1 / 34.4 / 34.6 | 19.9 / 16.2 / 16.8 | +3.8 / +5.3 / +3.6 % | −1.9 / −3.2 / −1.1 % |
| **T3 (seed 1 only)** | | | **+21 / +20 / +20 %** | −11 / −9 / −13 % |

**T3 is held. It breached the limit in the seed-1 smoke run, so the full set was not run and nothing was retuned.**
- Why it breaches: every variant keeps population-mean input fixed but redistributes it. The PV+ share of the mean input is 0.80 under T0, 1.00 under T1, 0.90 under T2 and 1.29 under T3. T3 moves input dorsally onto the PV+ cluster.
- T1 and T2 pass on 5 seeds. T1 is the closest at +8.4%.

## 1. Rate profiles by type (PD, graded depth; healthy in brackets)

| x | 0.05 | 0.25 | 0.45 | 0.65 | 0.85 | 0.95 |
|---|---|---|---|---|---|---|
| PV− T0 | 20.0 (24.9) | 18.5 | 17.1 | 15.2 | 13.2 | 12.3 (12.9) |
| PV− T1 | 15.6 (18.9) | 15.6 | 15.9 | **17.1 (21.0)** | 16.4 | 15.9 (18.9) |
| PV− T2 | 18.1 (22.2) | 17.3 | 16.5 | 15.7 | 13.9 | 12.9 (15.1) |
| PV+ T0 / T1 / T2 | 37 / 37 / 37 | | | 33 / 35 / 34 | | 30 / 36 / 33 |

- **T0.** The ventral→dorsal fall in PV− rate (about 1.6–1.9×) is the cortical gradient.
- **T2.** The fall is smaller (about 1.4×).
- **T1.** The cortical gradient is gone. The PV− profile is flat apart from a small hump at 0.65 (+1.5 Hz), where the PV+ cluster sits. I did not test the cause.
- **PV+ per-type profiles have n_eff < 8 at every bin except 0.65 (8.5).** The ventral bins hold almost no PV+ cells (n_eff 0–0.2), so the PV+ values there are kernel extrapolation. Treat PV+ profiles as descriptive only.
- **Uniform vs graded depth.** Rate profiles are almost the same at either depth. The pallidal depth changes beta, not mean rate.

## 2. Rate-normalised PD − healthy beta (×10⁻³)

| variant / depth | max (n_eff) | peak − mean vs SEM | seed maxima |
|---|---|---|---|
| T0 graded | 0.60 (21.5) | 6.5 vs 4.2 (1.5 SEM) | 0.75, 0.05, 0.65, 0.85, 0.95 |
| T0 uniform | **0.05 (7.9, flagged)** | 12.8 vs 11.0 (1.2) | 0.05, 0.85, 0.1, 0.05, 0.05 |
| T1 graded | 0.70 (20.4) | 10.1 vs 7.3 (1.4) | 0.4, 0.95, 0.85, 0.7, 0.7 |
| T1 uniform | 0.25 (19.5) | 5.0 vs 5.1 (1.0) | 0.05, 0.35, 0.95, 0.3, 0.9 |
| T2 graded | 0.65 (21.3) | 6.0 vs 7.1 (0.85) | 0.05, 0.65, 0.95, 0.3, 0.7 |
| T2 uniform | 0.70 (20.4) | 4.9 vs 8.5 (0.6) | 0.65, 0.3, 0.75, 0.05, 0.7 |

- **Graded depth.** The seed-mean maximum lands at 0.60–0.70 in every variant, where the graded pallidal drive is centred (0.65).
- **No maximum anywhere clears 2 SEM above the profile mean.**
- **Uniform depth.** There is no interior maximum.
  - Under T0 the maximum is at the ventral edge. It rests on n_eff 7.9 and the ventral cortical bias, so it is flagged.
  - Under T1 and T2 the profile is flat to within seed noise.
- **What the graded curves look like** (figure, row 2). The graded-depth curve rises from low values ventrally to a plateau or broad hump over 0.55–0.8. T0 also has a dorsal plateau, so the "peak" is a step more than a localised maximum.

## 3. PPC magnitude, PD − healthy

| | graded max (n_eff) | uniform max (n_eff) | mean PPC PV+ / PV− (graded; uniform) |
|---|---|---|---|
| T0 | 0.75 (18.7) | 0.05 (7.9, flagged) | 0.026 / 0.023; 0.042 / 0.074 |
| T1 | 0.70 (20.4) | 0.05 (7.9, flagged) | 0.020 / 0.030; 0.047 / 0.086 |
| T2 | 0.60 (21.5) | 0.25 (19.5) | 0.026 / 0.021; 0.042 / 0.070 |

- **PPC peaks agree with the normalised measure under graded depth** (within 0.05–0.15).
- **Uniform depth.** PPC is high and almost flat, decreasing dorsally. The 0.05 maxima rest on flagged bins.
- **PV− ≥ PV+ per spike.** This holds everywhere except under graded depth, where the two are about equal. That repeats the edge-artefact reanalysis.
- **T2 graded.** PPC falls below healthy at the dorsal edge (−0.036 at 0.95, n_eff 7.9, flagged).

## 4. Absolute modulated output ν·m (Hz)

| | PD − healthy max (n_eff) | peak − mean vs SEM | PD absolute max | PD ν·m at 0.05 / 0.65 / 0.95 |
|---|---|---|---|---|
| T0 graded | 0.60 (21.5) | 1.16 vs 0.48 (2.4 SEM) | 0.60 | 4.0 / 4.8 / 3.8 |
| T0 uniform | 0.60 (21.5) | 0.73 vs 0.58 (1.3) | **0.05 (flagged)** | 6.7 / 5.8 / 3.9 |
| T1 graded | 0.65 (21.3) | 0.97 vs 0.57 (1.7) | 0.65 | 2.9 / 6.4 / 5.1 |
| T1 uniform | 0.85 (13.6) | 0.36 vs 0.54 (0.7) | 0.65 | 5.6 / 7.8 / 6.7 |
| T2 graded | 0.65 (21.3) | 0.93 vs 0.74 (1.3) | 0.60 | 4.0 / 5.5 / 3.2 |
| T2 uniform | 0.65 (21.3) | 0.98 vs 0.61 (1.6) | 0.65 | 5.3 / 7.2 / 5.1 |

- **ν·m is the only readout with a maximum above 2 SEM:** T0 graded, PD − healthy at 0.60.
  - Its healthy ν·m falls dorsally (4.1 → 2.3 Hz) because the T0 cortical gradient lowers dorsal rates.
  - That makes the dorsal PD − healthy difference larger than the per-spike map alone.
- **Under T1 and T2 the absolute PD ν·m peaks at 0.60–0.65 at uniform depth too.** The per-spike map there is flat. So this absolute peak comes from rate, at the PV+ cluster and the PV− hump, not from beta localisation.
- **Uniform-depth absolute ν·m is higher than graded everywhere.** Uniform depth puts full pallidal beta into every cell.

## 5. Depth contrast: graded − uniform

The measure is [diff(0.65) − diff(0.25)] at graded depth minus the same at uniform depth, per seed, mean ± SE.

| variant | normalised (×10⁻³) | | ν·m (Hz) | |
|---|---|---|---|---|
| T0 | 17.6 ± 7.7 | 2.3 SE | 1.03 ± 0.67 | 1.5 SE |
| T1 | 26.0 ± 9.2 | 2.8 SE | 1.24 ± 0.68 | 1.8 SE |
| T2 | 3.6 ± 9.7 | 0.4 SE | −0.33 ± 0.78 | −0.4 SE |

- **T0 and T1: the per-spike depth contrast survives at 2.3–2.8 SE.**
- **T2: it does not survive (0.4 SE).**
  - Under T2, both depths give similar hump shapes in the normalised map at 0.65–0.70.
  - The graded seeds scatter: maxima at 0.05, 0.65, 0.95, 0.3 and 0.7.
  - I have not found a mechanism for this. It could be seed noise. With 5 seeds and SEs around 9 × 10⁻³, a true contrast as large as T0's (17.6) cannot be excluded.
- **ν·m contrast does not clear 2 SE in any variant.**
- **For reference** (normalised; `results/edge_artefacts/rate_normalised.md`):
  - configuration D original grid, graded vs uniform: 28.7 ± 7.8;
  - T0 with F1 + F2 here: 17.6 ± 7.7.

## Plain answer

1. **Removing the ventral cortical bias (T1) does not remove the graded-depth effect.** Per spike, the PD − healthy map still rises into 0.6–0.75, and the depth contrast over the uniform control survives (2.8 SE).
   - So the T0 rate gradient does not create that contrast.
   - The uniform-depth ventral-edge maximum under T0 does come from the cortical bias plus an edge bin (flagged). It disappears under T1 and T2.
2. **With the Jeon-fit conductance gradient (T2), the contrast is not resolved** (0.4 SE).
3. **No variant shows a localised per-spike maximum clearing 2 SEM.** What the graded drive produces is a rise into the dorsal half, centred near the drive's imposed centre (0.65).
   - That centre was set by hand in configuration D, so any effect here is the imposed drive showing through, not an emergent peak.
4. **The motor-territory variant (T3) could not be tested within the operating-point limit.**

## Jeon et al. 2022, Fig. 7H.2: what the Gria1 fraction is computed over

Source: Jeon et al., Cell Rep 38:110439 (2022), DOI 10.1016/j.celrep.2022.110439, PMID 35235786. I read the Methods and legends from the PDF. The publisher site was not reachable from this environment.

- **The fraction is over all detected cells, not within a cell type.**
  - Fig. 7H plots the "quantification of cells expressing each gene along STN axis 1", from n = 8,676 ± 315 cells in 6 STNs of 3 mice.
  - Cells are DAPI-segmented nuclei from RNAscope sections. Puncta are assigned to the nearest nucleus (≤ 1 µm). A support-vector machine, trained on 30 manually labelled cells per section, calls cells positive.
  - The same "total cells" denominator includes non-neuronal nuclei. Fig. 5 gives Pvalb+/Slc17a6+ as 10.2% of total cells.
- **Fig. 7I is the only within-type panel.** It shows the Gria1:Gabra1 ratio (lower in PV+ cells, higher in ventromedial PV− cells), not the Gria1-positive fraction along the axis.
- **Can PV+/PV− composition explain the gradient? Only partly.**
  - Along the axis the Gria1 fraction falls 0.55 → 0.28 while the Pvalb fraction rises about 0.07 → 0.3.
  - Even if every PV+ cell were Gria1-negative, composition would explain only about 0.55 → 0.41, roughly half of the drop.
  - The paper reports that most cells co-express Gria1 and Gabra1, so the real share is likely smaller.
- **A more likely confound, which is my inference and not tested in the paper:** the share of neurons among DAPI nuclei.
  - The authors note that the Gria1/Gabra1 gradients are "similar to the expression pattern of Slc17a6".
  - Both receptors decline dorsally together.
  - So the 7H.2 gradient may reflect neuronal density per nucleus rather than receptor expression per neuron.
  - That would argue against reading it as a per-cell conductance gradient, which is what T2 assumes.

## Caveats

- One wiring (loc 1, conn 4) and 5 input seeds; every seed-noise statement uses those 5 seeds.
- x = 0.05 and 0.95 rest on n_eff 7.9 (flagged). PV+ per-type profiles are below 8 everywhere except 0.65.
- Normalisation uses each cell's whole-run mean rate.
- The depth contrast uses two fixed positions (0.65, 0.25), chosen before these runs, as in the edge-artefact report.
