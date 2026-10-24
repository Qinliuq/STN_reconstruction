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

## Protection against depolarization block (2026-10-24): new defaults PV+ 1.2, PV- 1.15

Strong depolarizing steps (`--protocols depol`: +0.1 ... +0.8 nA for 4 s). The table shows rate,
spike trough and the end-of-step Vm at +0.8 nA. `*_depol.txt` holds all steps; `*_ref.txt` holds the
1.15 / 1.05 refinements; `*_fine.txt` holds the standard protocols at 1-1.3.

| Kv3.1 scale | PV+ at +0.8 nA | PV- at +0.8 nA |
|---|---|---|
| 0.5, 0.75 | **depolarization block** (silent at -27 mV after the first spikes) | **block** |
| 1.0 | 448 Hz small oscillation, trough -35 mV, AP width 0.17 ms (spikes do not repolarize: near-block) | 440 Hz, trough -35 mV (near-block) |
| 1.05 | | 429 Hz, trough -36 mV (near-block) |
| 1.1 | 425 Hz, trough -37 mV (near-block) | **219 Hz, full spikes, trough -56 mV** |
| 1.15 | **220 Hz, full spikes, trough -56 mV** | 218 Hz, trough -57 mV |
| 1.2 | 220 Hz, trough -57 mV | 217 Hz, trough -57 mV |
| 1.3 / 1.5 | 219 / 217 Hz, trough -57 / -59 mV | 216 / 215 Hz |

- Up to +0.5 nA every scale >= 0.75 fires full spikes (0.5 is already near-block at +0.5 nA).
- At +0.8 nA the transition is sharp. Below a threshold, spikes collapse into a small, non-repolarizing
  oscillation around -30 mV or a full block. Above it, the cell fires fully repolarized spikes at
  ~220 Hz.
- The threshold lies between 1.1 and 1.15 for PV+ and between 1.05 and 1.1 for PV-. The default of 1.0
  sits just below it in both types.
- The new defaults (`cfg.PVP_gKv31_scale = 1.2`, `cfg.PVN_gKv31_scale = 1.15`) are the bottom of the
  protective window plus a ~0.05 margin, keeping PV+ slightly above PV-.

Cost at the new defaults, from `*_fine.txt` (default 1.0 -> new):

| | PV+ (1.2) | PV- (1.15) |
|---|---|---|
| spontaneous rate | 17 -> 14 Hz | 11.3 -> ~9.7 Hz |
| +0.05 nA rate | 55 -> 53 Hz | 50 -> ~47.5 Hz |
| PV+ bursts at -0.16 nA | 9 x 5.3 -> 11 x 4.0 spikes | (none) |
| rebound | 14 -> 12 spikes, 161 Hz | 12 -> 7-9 spikes, 106-159 Hz |

PV+ keeps its rhythmic bursting and rebound; PV- keeps tonic firing. Both are slightly less excitable.
