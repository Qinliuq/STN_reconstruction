# Kv3.1 in isolated PV+ and PV- cells, 2026-10-18

Kv3.1 (`mod/Kv31.mod`, densities from `cells/sth-data/cell-gk_Kv31`) scaled per cell type with
`cfg.PVP_gKv31_scale` / `cfg.PVN_gKv31_scale` (default 1.0), swept with
`python src/single_cell.py --cell PVP --param gKv31 --values 0 0.25 0.5 0.75 1 1.5 2` (model state:
commit `2e7845f`; K_ATP off). `*.txt`: rate, ISI CV, bursting, spike width (time above -20 mV) and
voltage trough per protocol; rebound burst size and frequency.

## What Kv3.1 does

| Kv3.1 scale | PV+ spont. | PV+ spike width / trough | PV+ bursts at -0.16 nA | PV+ rebound | PV- spont. | PV- rebound |
|---|---|---|---|---|---|---|
| 0.25 | 32 Hz | 0.81 ms / -64.5 mV | 7 x 31 spikes | 2 spikes (depolarization block) | 22 Hz | 2 spikes |
| 0.5 | 27 Hz | 0.71 / -66.5 | 9 x 11 | 2 spikes | 18 Hz | 34 spikes, 302 Hz |
| 0.75 | 22 Hz | 0.67 / -68.0 | 9 x 8 | 22 spikes, 200 Hz | 14 Hz | 15, 182 Hz |
| **1 (default)** | **17 Hz** | **0.60 / -69.3** | **9 x 5.3** | **14, 164 Hz** | **11 Hz** | **12, 154 Hz** |
| 1.5 | 10 Hz | 0.56 / -71.4 | 11 x 3 | 8, 113 Hz | silent | 7, 81 Hz |
| 2 | 7 Hz | 0.50 / -73.0 | none | 8, 118 Hz | silent | 6, 105 Hz |

(Scale 0 is non-physiological: spikes fail to repolarize and cells go into depolarization block.)

- Kv3.1 speeds spike repolarization (narrower spikes, deeper trough) and lowers excitability:
  spontaneous rate falls steeply with Kv3.1, and PV- cells (weaker pacemakers) are silenced at 1.5x.
- High-frequency rebound bursts need an intermediate Kv3.1: too little and the rebound collapses into
  depolarization block, more than ~0.75x and rebound bursts get shorter and slower.
- PV+ rhythmic bursting under -0.16 nA shrinks with more Kv3.1 (31 -> 3 spikes/burst) and stops at 2x.

## Differentiating Kv3.1 between PV+ and PV- (PV+ slightly higher)

| setting | PV+ spont. | PV- spont. | PV+ bursts at -0.16 nA | PV+ rebound | PV- rebound |
|---|---|---|---|---|---|
| 1 / 1 (default) | 17 Hz | 11.3 Hz | 9 x 5.3 | 14, 164 Hz | 12, 154 Hz |
| PV+ 1.1 / PV- 0.9 | 15.3 Hz | 12.3 Hz | 10 x 4.5 | 14, 155 Hz | 12, 156 Hz |
| PV+ 1.2 / PV- 0.85 | 14.0 Hz | 12.7 Hz | 11 x 4.0 | 12, 161 Hz | 11, 167 Hz |

More Kv3.1 in PV+ than PV- makes PV+ less bursty and PV- slightly more excitable, so the two types
converge (rates, rebound size) instead of becoming more distinct. Kv3.1 differentiation therefore does
not make the burst behaviour more realistic; the scales stay at 1.0 and parkinsonian burstiness is
pursued through stronger GPe->STN synapses.
