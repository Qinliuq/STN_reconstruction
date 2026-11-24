# Data-derived PD drive topography: run specification as implemented (Priority 1, item 3)

Supersedes `spec_original.md` / `prompt_original.md` (kept for provenance) where the decisions of 2026-11-24 changed
them. Module: `src/stn_drive_topography.py` (the supplied module, with the changes listed in its docstring).
Depth scales: `params/drive_topography.json` (`tools/drive_topography/precompute.py`), committed before any run.
Sweeps: `sweeps/drive_topography_c.py` (variant c, controls, healthy references), `sweeps/drive_topography_ab.py`
(variants a, b). Analysis and verdict: `analyze.py` (criteria fixed before any run).

## Decisions (2026-11-24)

- **A. Seeds.** Input seeds 1–5 (`ctx_seed = s`, `gpe_seed = 100 + s`) on one fixed network (`seeds['loc'] = 1`,
  `seeds['conn'] = 4`), as in Table 3. Every peak-position result is conditional on this one wiring; no claim of wiring
  robustness. Follow-up (not blocking): conditions 1, 2, 3 with three extra `seeds['loc']` / `seeds['conn']` draws.
- **B. Beta depth.** Per GPe group from s(c_g) at its home position (configuration D: 12 even homes, σ 0.1), received by
  the cells through the existing projection; s renormalized so the connection-weighted mean depth received by the
  cells equals M. No new random draws (connectivity identical to configuration D).
- **C. Depth cap.** d_g = 0.1 + (M − 0.1)·s_g (0.1 = healthy depth), one M per comparison set, control included:
  M = 0.1 + 0.9 / max(s_g over the set's conditions and seeds), rounded down.

  | set | members | binding max s_g | M |
  |---|---|---|---|
  | advanced | 1, 2, 2b, 3, 4, 6 | 1.150 (2b) | **0.8825** |
  | early | 5, 5c (own uniform control) | 1.563 (5) | **0.6757** |

  Consequences: the advanced set runs at a ~12% weaker total beta drive than configuration D's uniform drive (depth
  1.0), so absolute beta is not comparable with Table 3 (the within-set comparison, which is the claim, is not
  affected). Condition 5 is on another scale: it is compared only with 5c, never directly with 2 / 3 / 6. Condition 1
  is therefore a new run (M 0.8825), not a reuse of `Dref_uni_pd` (depth 1.0), which stays a Table 3 reference.
- **D. The change is scaled, not the total.** GPe→STN synapse scale of cell i = 1 + (2.0 − 1)·s(x_i), s renormalized
  to mean 1 over the realised cells (asserted |mean − 1| < 1e-6 per run): advanced 1.70–2.14, early 1.42–2.54, never
  below healthy; mean exactly 2.0.
- **E. PV+ placement.** Condition 2 keeps configuration D's skew-normal; **2b** = advanced drive with the digitized
  *Pvalb* fraction (Jeon 2022 Fig. 7H) as placement. 6 = sigmoid ramp (`pv_shape_monotonic`) via
  `pvp_placement = 'profile'`. 12 / 58 PV+ throughout. `pv_density` is retired as a model input.
- **F. Healthy references.** One healthy set per placement (the healthy state does not depend on the drive):
  gradient = `Dref_healthy`, pvalb = `FU1_healthy` (both `sweeps/out/emergent_peak_followup`), uniform and monotonic
  new. Variant (c) first, then (a) and (b).

Realised PV+ positions (fixed network, 12 cells):

| placement | PV+ axis positions |
|---|---|
| gradient (skew-normal) | 0.48 0.51 0.52 0.52 0.53 0.57 0.59 0.66 0.67 0.73 0.74 0.89 |
| pvalb | 0.20 0.23 0.37 0.42 0.52 0.57 0.59 0.66 0.67 0.73 0.74 0.89 |
| uniform | 0.03 0.06 0.20 0.25 0.35 0.43 0.48 0.59 0.66 0.72 0.83 0.97 |
| monotonic | 0.36 0.39 0.53 0.58 0.61 0.65 0.68 0.74 0.76 0.81 0.82 0.97 |

## Conditions (PD; × 5 input seeds; variant c = both hooks, a = synapse only, b = depth only)

| # | label | drive | PV+ placement | set |
|---|---|---|---|---|
| 1 | `ctl_1_uniform_control_pd` | uniform (same in all variants) | gradient | advanced |
| 2 | `*_2_territory_adv_pd` | territory, advanced weights | gradient | advanced |
| 2b | `*_2b_territory_adv_pvalb_pd` | territory, advanced | measured *Pvalb* | advanced |
| 3 | `*_3_pv_uniform_pd` | territory, advanced | uniform | advanced |
| 4 | `*_4_mouse_map_pd` | mouse/Jeon intervals, advanced | gradient | advanced |
| 5 | `*_5_territory_early_pd` | territory, early weights | gradient | early |
| 5c | `ctl_5c_uniform_control_early_pd` | uniform | gradient | early |
| 6 | `*_6_pv_monotonic_pd` | territory, advanced | monotonic | advanced |

Runs: variant c 40 PD + 10 healthy; variants a + b 60 PD.

## Analytic illustration (not predictions about these runs)

The table "Predicted input-profile behaviour" in `spec_original.md` multiplies the drive by `pv_density` (Gaussian,
peak 0.68, σ 0.22), which the model never uses (its placements are listed above). It is an analytic illustration
of the drive × PV+ conjunction only; its peak positions are not predictions for these simulations, and its ratios are
not compared with simulated beta ratios.

## Readouts and acceptance criteria

As in `spec_original.md`, implemented in `analyze.py`. "Across 5 seeds" = the seed-mean PD − healthy difference map
(seed-matched healthy run, same placement); per-seed peaks reported next to it. The hypothesis survives only if all of
(1) condition 2 maximum in [0.55, 0.80], (2) condition 3 maximum > 0.85, (3) condition 6 maximum > 0.85, (4) condition
2's difference-map 0.65/0.25 ratio > condition 1's. 2b is reported side by side with 2 under rules 1 and 4.
