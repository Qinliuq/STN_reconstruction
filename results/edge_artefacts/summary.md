# Edge and interior maxima: artefact check — results, 2026-11-27

Configuration D, 58 cells, variant b (beta-depth hook only, M 0.8825), input seeds 1–5, one wiring (`loc` 1, `conn` 4).
Every PD condition is compared with a healthy reference that uses the same fixes and the same placement.
"Maximum" is the peak of the seed-mean PD − healthy map at positions with n_eff ≥ 5, as before; bins with n_eff < 8 are flagged.

- Step 1 diagnostics: `step1_diagnostics.md`.
- Runs: 150 new runs (R1: 5 conditions + geo; R2: 5 conditions; R3: 5 conditions + geo; R0 geo; matched healthy references). All succeeded. Nothing was retuned.
- Files: `edge_artefacts.csv` (one row per run × condition), `verdict_by_run.csv` (previous criteria per run), `difference_maps.png`.
- Deviation from the request: the geometry control was also run in R1, not only in R0 and R3.

## Operating point: one breach in R3 (healthy bursts), and R4 held

| run | PD rate PV+ / PV− (max abs change vs R0) | PD bursts | healthy rates | healthy bursts | within limits |
|---|---|---|---|---|---|
| R1 (F1) | ≤ 2.1% / ≤ 4.3% | −6.7 … +1.4% | ≤ 2.1% | +11.6 … +17.3% | yes |
| R2 (F2) | ≤ 1.3% / ≤ 1.4% | +4.7 … +16.0% | ≤ 1.1% | +9.6 … +20.5% | yes |
| R3 (F1+F2) | ≤ 1.5% / ≤ 3.1% | +6.8 … +14.9% | ≤ 1.6% | **+20.0 … +33.7%** | **no** for the gradient, monotonic and *Pvalb* placements (conditions 1, 2, 6, 2b) |
| R4 (F1+F2+F3) | seed-1 smoke run: **−17% / −30%** | **−53%** | — | — | **no**: held, not run |

- **R3.** The breach is in the healthy reference: its burst fraction rises from 0.046–0.049 to 0.059–0.063. All PD rates and PD bursts stay within limits. As instructed, I did not recalibrate.
  - The R3 maps below are shown for completeness only.
  - Any R3 result in conditions 1, 2, 6 and 2b compares PD with a healthy state that has moved, so none of them is used as evidence. The uniform-placement R3 runs (condition 3, geo) are within limits.
- **R4.** In the seed-1 smoke run, delivered GPe input was 147 vs 165 nS/s and total weight per cell was 5.27 vs 5.45 nS.
  - With 48 groups, each cell receives ~12 groups instead of ~3. Inhibition is smoother, and rebound bursting collapses.
- **Delivered GPe input** is logged for every new run: 160–165 nS/s population mean in every run and condition, within 3%. Ventral and dorsal means are in the CSV.
  - The R0 runs reused from earlier sweeps predate the logging; only R0 geo has it.

## 1. Does the 0.05 maximum survive F1 (centred grid, reflected edges)?

**Not for the non-uniform PV+ placements; yes for the uniform ones.**

| condition | R0 max | R1 (F1) max | note |
|---|---|---|---|
| 1 uniform control | 0.05 | **0.70** | edge gone; interior hump moved from 0.55 to 0.70 |
| 2 territory, skew-normal | 0.60 | 0.65 | edge gone (6.0 vs 12.4 interior) |
| 2b territory, *Pvalb* | 0.05 | 0.45 | edge gone; flat map (range 3.6–7.0) |
| 6 monotonic PV+ | 0.60 | 0.55 | edge 1.4 (was 7.6) |
| 3 uniform PV+ | 0.05 | **0.05** | survives (n_eff 8.1) |
| geo (uniform depth + PV+) | 0.05 | **0.05** | survives (17.1 vs 3.3 dorsal; n_eff 8.1) |

- **Where it disappears.** With the gradient, monotonic and *Pvalb* placements, the ventral-edge maximum disappears once positions stop being clipped and the grid is symmetric.
  - For these placements it was a property of the asymmetric, clipped grid together with the hump-shaped healthy reference (Step 1).
- **Where it survives.** With PV+ cells placed uniformly (condition 3, geo), the edge stays the maximum under F1. In R3 geo (within limits) it is larger still: 28.2 vs 5.7.
  - Its n_eff is 8.1, just above the flag line.
  - It is not a clipping or coverage artefact. I have no explanation for it from these runs.
- **Previous criteria under F1.** No criterion changes its pass/fail result: criterion 1 pass, criteria 2 and 3 fail, criterion 4 pass (2.03 vs 1.71). Rules 1 and 4 for 2b still fail. **The hypothesis still fails.**
- **What it undermines.** Criterion 1 does not discriminate under F1. The uniform-drive control (condition 1) now peaks at 0.70, inside the [0.55, 0.80] window, with no territory drive.
  - The earlier "condition 2 peaks at 0.60" result is therefore no evidence for the drive. The interior hump in conditions 1 and 2 follows the PV+ cluster, which F1 shifts by +1/24 (gradient now 0.52–0.78).
  - Condition 1's interior hump moved 0.55 → 0.70 with it. This matches Step 1: PD shared beta tracks rate², and PV+ cells fire faster.
  - Condition 1's 0.65/0.25 ratio rises from 1.27 to 1.71, so the bar for criterion 4 is now higher.

## 2. Does the ~0.6 interior maximum survive F2 (per-cell GPe weight normalised)?

**Yes, shifted 0.05 ventrally. It is not a GPe-coverage effect, but it becomes less condition-specific.**

| condition | R0 max | R2 (F2) max | 0.65/0.25 ratio R0 → R2 |
|---|---|---|---|
| 1 uniform control | 0.05 (interior 0.55) | 0.55 | 1.27 → 1.52 |
| 2 territory, skew-normal | 0.60 | 0.55 | 2.11 → 1.90 |
| 6 monotonic PV+ | 0.60 | 0.55 | 1.71 → 1.39 |
| 3 uniform PV+ | 0.05 (interior 0.30) | 0.50 (shallow: 11.4 → 13.5) | 0.84 → 1.02 |
| 2b territory, *Pvalb* | 0.05 | 0.15 (11.4; second peak 0.55: 11.3) | 0.87 → 0.98 |

- Removing per-cell differences in total GPe weight leaves a hump at 0.50–0.55 in every condition, including uniform PV+ (condition 3), where R0 had none.
- Condition 2's contrast over the control drops from 1.66× to 1.25× the control's ratio. Criterion 4 still passes.
- So the interior maximum is not created by uneven GPe coverage.
- Part of what separated condition 2 from condition 3 was coverage. Once coverage is equal, the hump sits near 0.5–0.55 regardless of PV+ placement or drive.
- With F2, every condition's dorsal end drops close to zero difference (condition 2: 0.4 at 0.95; flagged n_eff < 8).
- R3 (F1+F2) also puts the maximum at 0.55 in conditions 1 and 2, but those runs breach the healthy-burst limit (see above) and are not used.

## 3. Does F3 (48 GPe groups) move the peak, or only the jitter?

**Not answered.** R4 breaks the operating point (PV− rate −30%, bursts −53% in the smoke run). Under the no-retune constraint it was held rather than recalibrated.

To answer it, a 48-group drive would need to keep the ~3-groups-per-cell arrival statistics, for example by narrowing σ to keep the per-cell group count. That changes the drive statistics, so it is your call.

## 4. Does PPC tell the same spatial story as cross-power?

**No.** Locking to the imposed GPe beta (across-episode PPC₁) is near-universal in PD:
- 81–98% of cells are significantly locked;
- 92–100% of PV+ cells and 76–98% of PV− cells;
- healthy: 0–12%, around the 5% false-positive rate.

Its PD − healthy map varies little (≈ 0.02–0.09) and its maximum rarely coincides with cross-power's:

| run / condition | cross-power max | PPC max | PPC n_eff |
|---|---|---|---|
| R0 / 1 | 0.05 | 0.95 | **5.7 (< 8)** |
| R0 / 2 | 0.60 | 0.95 | **5.7 (< 8)** |
| R0 / 3 | 0.05 | 0.85 | 11.2 |
| R0 / 6 | 0.60 | 0.60 | 21.0 |
| R0 / 2b | 0.05 | 0.95 | **5.7 (< 8)** |
| R1 / 1 | 0.70 | 0.20 | 17.5 |
| R1 / 2 | 0.65 | 0.95 | **7.9 (< 8)** |
| R1 / 3 | 0.05 | 0.55 | 21.4 |
| R1 / 6 | 0.55 | 0.85 | 14.0 |
| R1 / 2b | 0.45 | 0.55 | 21.4 |
| R2 / 1 | 0.55 | 0.35 | 21.0 |
| R2 / 2 | 0.55 | 0.40 | 21.0 |
| R2 / 3 | 0.50 | 0.45 | 21.6 |
| R2 / 6 | 0.55 | 0.55 | 21.1 |
| R2 / 2b | 0.15 | 0.75 | 16.7 |

- The maxima agree within 0.05 in 4 of 15 cases (R0/6, R2/3, R2/6, and R1/2b at 0.1 apart).
- Under R0, PPC rises toward the dorsal end in conditions 1, 2 and 2b, but those maxima sit on bins with n_eff 5.7. **Any "dorsal locking" reading of R0 hinges on n_eff < 8 and is flagged.**
- With F2, the PPC and cross-power maps line up more closely in shape (R2 row of the figure), though the PPC maxima still sit 0.15–0.20 more ventral in conditions 1 and 2.

**Interpretation.** Pairwise cross-power is rate²-weighted, so its spatial structure is dominated by where fast-firing (PV+) cells sit. Rate-independent locking says nearly every cell follows the imposed beta and shows no stable spatial maximum.

## What this changes relative to earlier results

1. **The ventral-edge maxima are grid artefacts** for the gradient, monotonic and *Pvalb* placements. This was the frequent "0.05" in the drive-topography tables, including condition 2b and the uniform controls.
   - The edge in uniform-PV+ placements is not an artefact and remains unexplained.
2. **The interior PD beta peak is a PV+-cluster / rate effect, not a drive effect.**
   - Under F1 the uniform-drive control peaks inside the hypothesis window (0.70).
   - Under F2 a hump at 0.50–0.55 appears in every condition.
   - This weakens the drive-topography "criterion 1 pass (0.60)" in variant b, and earlier claims that a peak at ~0.6–0.7 reflects the drive.
3. **The negative verdict on the drive-topography hypothesis stands** under F1 and F2: criteria 2 and 3 fail in every run. No run puts the maximum dorsal of 0.85.
4. **Cross-power is the wrong readout for spatial beta in this model** if the question is locking rather than rate-weighted co-fluctuation. The two disagree on position in most cases.

## Caveats

- All results use one wiring (`loc` 1, `conn` 4) and 5 input seeds. PPC is pooled over seeds per cell.
- R3 (conditions 1, 2, 6, 2b) breaches the healthy-burst limit and is not used for conclusions. R4 was not run.
- The F1 edge bins have n_eff ≈ 8.1 in uniform placements, just above the flag line.
- The 0.95 bins are below 8 under the original grid (5.7) and in some F1 maps.
