# Run spec — data-derived drive topography (Priority 1, item 3)

Companion to `stn_drive_topography.py` and to the project doc
`claude/stn-denervation-territory-weights.md`.

Starting point: configuration D (58 STN cells, equal Kir, PD GPe→STN scale 2.0,
PV+/PV− GABA:glutamate ratio). Axis x ∈ [0,1], 0 = anteromedial/ventromedial,
1 = posterolateral/dorsolateral.

## Wiring it in

Two hook points, both per cell, from the cell's normalised position `x` on the
principal axis:

```python
from stn_drive_topography import drive_scaling, pv_density, CONDITIONS, SEEDS

s = drive_scaling(x, **CONDITIONS[cond]["drive"])    # mean 1 over the population
gpe_scale_cell  = 2.0 * s        # was uniform 2.0
beta_depth_cell = 1.0 * s        # was uniform 1.0
```

Run each of (a) GPe→STN potentiation only, (b) beta depth only, (c) both —
the three variants of step 3 in the doc. Start with (c).

`drive_scaling` is mean-normalised, so every condition delivers the same **total**
PD change as the uniform control; only its distribution over x differs. Keep the
healthy state exactly as measured, including its ventral cortical bias — the
claim is about the PD − healthy contrast, not about absolute beta.

## The six conditions (× 5 connectivity seeds = 30 runs per variant)

| # | Condition | Drive | PV+ placement | Tests |
|---|---|---|---|---|
| 1 | `1_uniform_control` | uniform, same total | gradient | control — reproduces the "uniform drive" cell of the 2×2 |
| 2 | `2_territory_adv` | territory-weighted, advanced weights | gradient | the hypothesis |
| 3 | `3_pv_uniform` | territory-weighted, advanced | uniform | is the peak position set by PV+? |
| 4 | `4_mouse_map` | territory-weighted, mouse/Jeon map | gradient | species sensitivity of the territory map |
| 5 | `5_territory_early` | territory-weighted, early weights | gradient | does disease stage change contrast or location? |
| 6 | `6_pv_monotonic` | territory-weighted, advanced | monotonic PV+ | is the interior peak due to PV+ falloff? |

## Predicted input-profile behaviour (before any simulation)

From `stn_drive_topography.py` directly — these are **input** profiles, not beta:

| # | drive peak | product peak | product ratio 0.65/0.25 |
|---|---|---|---|
| 1 | — (flat) | 0.65 | 6.69 |
| 2 | 0.95 | 0.65 | 8.22 |
| 3 | 0.95 | **0.95** | 1.23 |
| 4 | 0.95 | 0.65 | 8.26 |
| 5 | 0.95 | 0.75 | 14.04 |
| 6 | 0.95 | **0.95** | 8.63 |

Read these for **peak position**. The network compresses contrast heavily — in the
2×2 a far steeper input gradient yielded a beta ratio of only 2.43 — so do not
compare these ratios against beta ratios.

## Readouts (per run, then pooled over 5 seeds)

- PD − healthy beta **difference** map along x; peak position per seed
- 0.65/0.25 ratio on both the difference map and absolute beta
- STN beta power, synchrony, fraction of spikes in bursts (as in Table 3)
- pairwise shared beta by cell type (as in Table 5) — does the PV+ : PV− ratio
  of ~2.4–3.4 survive a data-derived drive?

## What each comparison decides

- **2 vs 1** — does a data-derived topography localise beta at all, given the
  same total PD change? If 2 ≈ 1 the whole territory-weighting idea is dead.
- **2 vs 3** — expected to be the decisive pair. Prediction: the peak moves to
  the dorsolateral edge without the PV+ gradient.
- **2 vs 6** — isolates the *falloff* of PV+ from its mere presence. If 6 also
  peaks at 0.65 the explanation is not the falloff and the account needs redoing.
- **2 vs 4** — the input profiles barely differ (peak 0.65 either way), so a
  difference in simulated beta would mean the network, not the map, is sensitive
  to the species assumption.
- **2 vs 5** — prediction: contrast grows, peak does not move. If the peak moves,
  position tracks disease stage, which is testable against staged patient LFPs.

## Acceptance criteria

The hypothesis survives if, across 5 seeds:
1. condition 2 puts the PD − healthy beta maximum at x = 0.55–0.80, and
2. condition 3 moves it dorsolateral of 0.85, and
3. condition 6 moves it dorsolateral of 0.85, and
4. condition 2's difference-map ratio exceeds condition 1's.

If 1 holds but 3 and 6 do not, the peak is being set by something else in the
network (geometry, boundary effects, connectivity density) and that needs
identifying before the result is reported.

## Known soft spots

- `PV_PEAK` and `PV_SIGMA` are assumptions — Jeon gives no percentages, only
  "dorsolateral and middle". Sweep peak 0.55–0.80 and sigma 0.15–0.30; the
  conjunction result depends on them. (To-do item 8: Allen MERFISH could set
  these properly.)
- The limbic weight (0.45 advanced / 0.20 early) has no direct source; sweep
  0.3–0.6 and 0.1–0.3.
- `n` and `K` in the saturating transfer are unconstrained; the doc sweeps
  n ∈ {1,2,4}, K ∈ {0.5,0.7,0.85}. Peak position was insensitive to these;
  contrast was not.
- Territory map is human (Plantinga), PV+ gradient is mouse (Jeon), recordings
  are mouse. State this in the Methods rather than letting a reviewer find it.
