# Kir in isolated PV+ and PV- cells: removing the 5x asymmetry, 2026-10-25

Previously PV+ had Kir scale 0.5 and PV- 0.1, so a uniform Kir block perturbed PV+ ~5x more. Here
both types are swept over the same range with
`python src/single_cell.py --cell PVP --param gKir --values 0.1 0.2 0.3 0.4 0.5`, and likewise for PVN.
Other scales are as in `src/cfg.py` (Kv3.1 PV+ 1.2 / PV- 1.15). `*.txt` hold the full protocol
output; the figures show the voltage traces.

| Kir scale | PV+ spont. | PV+ at -0.16 nA | PV+ rebound | PV- spont. | PV- at -0.16 nA | PV- rebound |
|---|---|---|---|---|---|---|
| 0.1 | 32 Hz tonic | 14 bursts x 3 spikes | 12 spikes, 154 Hz | 10 Hz tonic | silent | 8 spikes, 156 Hz |
| 0.2 | 27 Hz | 13 x 3 | 12, 156 Hz | 7.7 Hz | silent | 8, 159 Hz |
| **0.3** | **22 Hz** | **12 x 3** | **12, 159 Hz** | **6.7 Hz** | **silent** | **8, 160 Hz** |
| 0.4 | 18 Hz | 11 x 3 | 12, 161 Hz | silent | silent | 8, 160 Hz |
| 0.5 | 14 Hz | 11 x 4 | 12, 161 Hz | ~1 Hz | silent | 8, 165 Hz |

- The bursting phenotype does not depend on Kir.
  - PV+ bursts rhythmically under -0.16 nA and fires a 12-spike rebound at every Kir level.
  - PV- never bursts and fires an 8-spike rebound at every level.
  - The phenotype difference is carried by CaT, CaL, HCN and SK, not by Kir.
- Kir mainly sets the spontaneous pacemaking rate. Each +0.1 lowers it by ~4-5 Hz in PV+. PV-, the
  weaker pacemaker, falls silent from 0.4.
- **New default: Kir 0.3 in both types** (`cfg.PVP_gKir_scale = cfg.PVN_gKir_scale = 0.3`), the mean
  of the earlier 0.5 / 0.1, so total STN Kir is about unchanged.
  - Both phenotypes are kept: PV+ is a 22 Hz pacemaker that bursts under hyperpolarization, and PV-
    is a 6.7 Hz tonic, non-bursting pacemaker.
  - PV+ is now faster at rest than before (14 Hz), and PV- slower (10 Hz).
