# Rate-normalised reanalysis of the stored spikes, 2026-11-28

No new simulations and no parameter changes. Script: `rate_normalised.py`.

Data:
- the R0–R3 runs of the edge-artefact set (variant b, 5 input seeds, matched healthy references);
- configuration D runs `Dref_healthy`, `Dref_pd` and `Dref_uni_pd` (`sweeps/out/emergent_peak_followup`), which are Table 5's three rows.

Outputs:
- `rate_normalised.csv`: maxima, ratios, rate and edge diagnostics per run × condition;
- `rate_normalised_pertype.csv`: Table 5-style per-pair values;
- `rate_normalised_profiles.json`;
- `rate_normalised.png`: rows R0–R3 show rate-normalised PD − healthy (black) and PPC magnitude PD − healthy (blue); the bottom row shows R0 rate profiles by type.

**Measures.**
- **Rate-normalised pairwise beta.** Each cell's 1 ms spike-count train is divided by its own mean count, so the pairwise 13–30 Hz cross-power becomes C_ij / (ν_i ν_j). This is the product of the two cells' shared fractional beta modulation: modulation per spike, independent of rate. The estimator is unchanged: weighted mean minus auto terms, kernel 0.12, n_eff in cells.
  - Raw cross-power recomputed this way reproduces Table 5: PD graded 17.2 / 5.1 / 9.7 ×10⁻⁶ for PV+ / PV− / cross pairs.
- **PPC.** Across-episode PPC₁ magnitude, pooled per cell over the 5 seeds; no fraction-locked readout is used.
- **Seed-noise check.** For every maximum I compare the peak height above the profile mean with the standard error across seeds at the peak (SEM).

**Flags.**
- R3 conditions 1, 2, 6 and 2b breach the healthy-burst limit (previous report). They are shown but not used for conclusions.
- R2 has no geometry control.

## 1. Firing-rate profiles by type

In every condition and state, the rate of each type falls from ventral to dorsal:

| | ventral end (0.05) | dorsal end (0.95) |
|---|---|---|
| PV− | ~20 Hz | ~11 Hz |
| PV+ | 35–44 Hz | ~29 Hz |

- This is the cortical input gradient (`ctx_topography` 0.65 − 0.5·x; a rate scale of 1.45× at the ventral edge). It is the same in every run.
- PV+ cells fire about 2× faster than PV− cells at every position. So ν² differs about 4× between the types, and about 3.6× along the axis within PV−.
- With the clustered PV+ placements (gradient, *Pvalb*), the whole-population rate profile has a hump at the PV+ cluster (~0.55). The raw PD cross-power profile tracks the kernel mean of ν² (r = 0.55–0.99 across runs).
- **The PV+ rate profile rests on n_eff < 8 at every position.** There are only 12 PV+ cells: per-type n_eff is 0–8.5 for the gradient placement and 1.8–3.8 for the uniform one. The PV− profiles (n_eff 8–21) are not flagged.

## 2. Table 5 and the 2×2, rate-normalised (×10⁻³ per pair; seed mean)

**Table 5 (configuration D)**

| row | PV+ pairs | PV− pairs | PV+–PV− | PV+ / PV− normalised | PV+ / PV− raw |
|---|---|---|---|---|---|
| healthy | 5.4 ± 1.3 | 6.9 ± 5.4 | 4.4 | 0.78 | 2.61 |
| PD, graded drive | 15.8 ± 8.6 | 21.1 ± 3.9 | 18.9 | **0.75** | 3.35 |
| PD, uniform drive | 32.2 ± 13.9 | 54.1 ± 20.3 | 39.8 | **0.59** | 2.40 |

**The original 2×2 cannot be rebuilt in full.** Only two of its cells exist on disk: graded × PV+ gradient (`Dref_pd`) and uniform × PV+ gradient (`Dref_uni_pd`), both above. The graded × uniform-PV+ and uniform × uniform-PV+ runs were never stored and were not re-simulated.

**Variant-b 2×2 (complete in R0 and R1)**

Each cell gives PV+ / PV− normalised, then raw.

| | PV+ gradient | PV+ uniform |
|---|---|---|
| territory drive, R0 | 22.7 / 35.3 = 0.64 (raw 2.69) | 21.3 / 58.3 = 0.36 (raw 1.62) |
| uniform drive, R0 | 18.4 / 40.9 = 0.45 (raw 1.95) | 19.2 / 44.1 = 0.43 (raw 1.97) |
| territory drive, R1 | 28.8 / 34.2 = 0.84 (raw 3.18) | 19.6 / 53.7 = 0.36 (raw 1.81) |
| uniform drive, R1 | 22.1 / 36.8 = 0.60 (raw 2.41) | 18.0 / 42.4 = 0.42 (raw 2.01) |

**Are PV+ pairs more strongly modulated per spike than PV− pairs? No, the opposite.**
- PV− pairs carry more shared beta per spike in every PD condition: seed-mean ratio 0.34–0.84. Per seed, PV+ < PV− in 119 / 125 PD seed-runs. Healthy: 0.78, but only 3 / 5 seeds and PV− SD ±5.4, so that comparison is undecided.
- **PPC agrees** (per cell, rate-independent): mean PV− PPC exceeds PV+ in every condition (PV− 0.052–0.083 vs PV+ 0.025–0.050). In configuration D's graded drive the two are nearly equal (0.027 vs 0.025).
- **The raw 1.6–3.4× PV+ advantage is a rate effect.** The ν² ratio is about 4×. Once it is divided out, PV+ cells share less beta per spike.
- **PD gain per spike** (PD / healthy) is similar for the two types in Table 5: 2.96 vs 3.06 (graded), 6.0 vs 7.9 (uniform). In variant b it varies with placement (PV+ 3.4–22, PV− 3.1–10). Those ratios divide by small and noisy healthy PV+ values.

**Caveat on the interpretation.** At ~33 Hz, PV+ cells fire more than once per beta cycle, which caps their fractional modulation. "Less modulated per spike" is a statement about the readout, not a mechanism. Even so, both per-spike measures (normalised cross-power and PPC) agree.

## 3. PPC magnitude profile (PD − healthy)

PPC maxima per condition are in `rate_normalised.csv` and on the blue curves in the figure. Where both the rate-normalised and PPC profiles have a clear maximum, they now agree more often than raw cross-power and PPC did:
- R0 condition 6: 0.60 / 0.60;
- R1 condition 3: 0.50 / 0.55;
- R2 condition 6: 0.55 / 0.55;
- configuration D graded drive: **0.75 / 0.75**;
- R1 condition 2: 0.75 / 0.95 (n_eff 7.9, **flagged**).

Elsewhere they differ by 0.15–0.6.

**The original-grid PPC maxima at 0.95 (R0 conditions 1, 2, 2b) rest on n_eff 5.7 and are flagged.** So is the R0 condition 2b normalised maximum (0.95, n_eff 5.7).

## 4. The uniform-PV+ maximum at 0.05: explained by rate

- **Raw cross-power.** In the uniform-PV+ placements (condition 3, geo), 2 of the 5–7 cells within 0.15 of the ventral edge are PV+ (26–29%), against 0% in the gradient placements.
  - These ventral PV+ cells are the fastest cells in the network: 42–44 Hz, because of the cortical gradient.
  - So the edge has the highest kernel-mean ν² of any position. The PD population rate there is 26–27 Hz, against ~19 Hz at mid-axis and ~20 Hz at the same edge with the gradient placement.
  - The raw PD map follows ν² (r = 0.71–0.99).
- **Rate-normalised.** The edge maximum disappears in the uniform-PV+ condition 3 in every run (R0 0.65, R1 0.50, R2 0.60, R3 0.70) and in R0 geo (0.35). Under normalisation the edge is at or below the profile mean.
- **Two F1-grid exceptions are not cleanly removed.**

  | run | normalised max | n_eff | peak − mean vs SEM | seed maxima |
  |---|---|---|---|---|
  | R1 geo | 0.05 | 8.1 (borderline) | only 0.5 × 10⁻³ above the 0.45 hump | 0.05, 0.95, 0.4, 0.05, 0.3 |
  | R3 geo | 0.10 | 10.9 | 11.7 vs SEM 9.2 | 0.05, 0.1, 0.3, 0.05, 0.8 |

  R3 geo is within limits. In both, a weak ventral per-spike excess remains, comparable to seed noise.
- **Answer.** Yes, the uniform-PV+ 0.05 maximum is explained by the rate profile. The ventral cortical bias raises the rate of both types, and the uniform placement puts the fastest (PV+) cells at the ventral edge. In the uniform-drive, uniform-PV+ F1 runs, a residual ventral excess per spike is at the noise level.

## 5. Plain statement

**Q1. Under the rate-normalised measure, does any condition show a PD − healthy spatial maximum?**

**Only the original configuration-D graded drive does.**
- **Configuration D graded drive** (`Dref_pd`: depth 1.0 Gaussian at 0.65, width 0.15):
  - normalised maximum at 0.75 (n_eff 16.8);
  - peak − mean 12.5 ×10⁻³ against seed SEM 5.4 (2.3 SEM);
  - 3 / 5 seeds peak at 0.65–0.75;
  - PPC peaks at the same place (0.75).
- **Variant-b conditions inside the limits (R0–R2, R3 conditions 3 and geo):**
  - the mean-profile maximum is at most 1.65 SEM above the profile mean, and usually below 1.1;
  - the one exception is R1 geo, at 1.9 SEM. It sits at the ventral edge (n_eff 8.1, borderline) and is essentially tied with a hump at 0.45;
  - per-seed maxima scatter across the axis, e.g. R0 condition 2: 0.6, 0.8, 0.95, 0.7, 0.7; R2 condition 2: 0.6, 0.65, 0.5, 0.1, 0.2.
  - The seed-mean maxima land at 0.55–0.75 for the territory drive, in both PV+ placements. But none clears seed noise as a localised maximum.
  - The uniform-drive control (condition 1) has no consistent maximum: 0.05, 0.75 and 0.35 in R0–R2.

**Q2. Does the depth-driven contrast over the uniform control survive?**

Measure: [diff(0.65) − diff(0.25)] for condition 2 minus the same for condition 1, per seed (×10⁻³, mean ± SE).

| run | contrast | about | survives? |
|---|---|---|---|
| R0 (original grid) | 12.0 ± 5.5 | 2.2 SE | yes |
| R1 (F1, centred grid) | 2.6 ± 8.8 | 0.3 SE | no |
| R2 (F2, coverage-normalised) | 8.1 ± 7.3 | 1.1 SE | no |
| configuration D, graded vs uniform drive | 28.7 ± 7.8 | 3.7 SE | yes |

- **For the variant-b (territory) drive it survives only on the original grid, at about 2 SE. It does not survive either bias fix.**
- **The configuration-D graded drive's contrast survives clearly.** That drive was hand-set at 0.65, so this is the imposed drive showing through, not an emergent effect.
- Per spike, the territory drive's contrast does not depend on the PV+ cluster: condition 3 (uniform PV+) gives an R0 contrast similar to condition 2 (per-seed mean 9.7 vs 13.0).

**What this undermines.**
1. **"PV+ cells are intrinsic beta amplifiers" / "PV+ pairs share 2.4–3.4× the beta" (`results/gaba_ratio`, Table 5).** True for raw cross-power and false per spike. The PV+ excess is their higher rate (ν² ≈ 4×). Both per-spike measures put PV− at or above PV+.
2. **"The upper-middle peak needs the PV+ gradient" (Table 3).** Under the rate-normalised measure, the graded drive's peak does not need PV+ cells. The raw-measure dependence on PV+ placement is largely where the fast cells sit.
   - Only the graded × gradient cell of that 2×2 was re-checkable: the PV+-uniform cells are not on disk.
3. **Drive topography, variant b.** The condition 2 interior peak and its contrast over the control are within seed noise per spike once either bias fix is applied.

## Caveats

- One wiring (`loc` 1, `conn` 4), with 5 input seeds per condition. The seed-noise statements use those 5 seeds.
- The normalisation divides by each cell's own mean rate over the whole 3 s, not the in-episode rate.
- PV+ per-type rate profiles and some 0.95 bins rest on n_eff < 8 (flagged above). No conclusion above depends on them except the edge-rate figures for PV+ cells, which are per-type means of 2–3 cells.
