# Kv3.1 x K_ATP in isolated PV+ cells: diagnostic and admissible region, 2026-11-18

This is a side test of sufficiency. Configuration D is unchanged: PV+ Kv3.1 1.2, PV- 1.15 and
Kir 0.3 / 0.3 in all states, and K_ATP is not part of it.

The cells are isolated PV+ cells with k_Ca = 20 uM, pdscale 1 and the other conductances as in
`src/cfg.py`. The model state is commit `0bef021`. Each cell is run for 11 s through the
standard protocols (statistics over 1-11 s), plus the 4 s depolarizing set. The admissibility
criteria are in `criteria.md` and were committed before the runs.

```
python src/single_cell.py --cell PVP --param katp_gbar --values 0 1e-5 2e-5 3e-5 5e-5 --duration 11000 \
    --save raw/kv<K>_kir<R>.npz --cfg katp_kca=0.02 PVP_gKv31_scale=<K> PVP_gKir_scale=<R>
python src/single_cell.py ... --protocols depol --save raw/kv<K>_kir<R>_depol.npz
python results/katp_kv31_map/analyze.py
```

Files:
- `map.csv`: the Kir 0.3 grid. `diagnostic.csv`: the Kv3.1 1.0 / 1.2 x Kir 0.3 / 0.5 cells.
- `map.png`: the map, with flags per cell.
- `analyze_output.txt`: the full table.
- `raw/`: spike times, logs and voltage traces per Kv3.1 / Kir setting.

## 1. Diagnostic: which change removed the K_ATP bursting?

Spontaneous firing at K_ATP 3e-5 S/cm2, in the 1-4 s window of the 2026-10-18 run:

| PV+ Kv3.1 / Kir | rate | ISI CV | bursts (>= 3 spikes) | +0.8 nA, 4 s |
|---|---|---|---|---|
| **1.0 / 0.5** (2026-10-18 cell) | 16 Hz | **1.14** | **5.3/s x 3.0 spikes** (12-14 ms ISIs, ~185 ms apart) | near-block (trough -36 mV) |
| 1.0 / 0.3 (Kir lowered only) | 12 Hz | 0.77 | none: doublets, 20 ms apart, every ~165 ms | near-block |
| 1.2 / 0.5 (Kv3.1 raised only) | 10 Hz | 0.88 | none (0.3/s at 5e-5) | full spikes (-58 mV) |
| 1.2 / 0.3 (configuration D) | 5.7 Hz | 0.01 | none: tonic | full spikes |

- **The 2026-10-18 bursting is reproduced**: 5.3 bursts/s of 3 spikes, CV 1.14, against ~5 Hz x 3
  and CV 1.1 before. The same cell bursts at 5e-5 too (4.3/s x 3, CV 1.25).
- **Each of the two configuration-D changes removes the bursting on its own.**
  - Lowering Kir 0.5 -> 0.3 at Kv3.1 1.0 turns the triplets into doublets (CV 0.77).
  - Raising Kv3.1 1.0 -> 1.2 at Kir 0.5 leaves at most rare triplets (0.3/s at 5e-5).
  - It is not a Kv3.1 effect alone. **Kir 0.5 is required**: at Kir 0.3, no Kv3.1 down to 0.9 gives
    regular bursting; the most is 0.2 bursts/s at Kv3.1 0.9 and 5e-5.
- **The bursting cell fails depolarization-block protection.** At Kv3.1 1.0 its spikes no longer
  repolarize at +0.8 nA (trough -36 mV), consistent with `results/kv31_single_cell/`.

## 2. Map (Kir 0.3): no admissible cell

| criterion | cells passing (of 30) |
|---|---|
| 1. spontaneous burst/pause | **0** (Poisson-surprise burst share 0 in every cell) |
| 2. no block at +0.8 nA | 14: Kv3.1 >= 1.15 at every gbar; Kv3.1 1.1 only with gbar >= 2e-5 |
| 3. -0.16 nA bursting + rebound kept | 18: every Kv3.1 at gbar <= 2e-5. From 3e-5 the -0.16 nA bursting is lost, while the rebound (10-15 spikes, > 150 Hz) survives everywhere |
| all three | **0** |

- **What K_ATP does at Kir 0.3.** It slows pacemaking (22-29 Hz -> 5-14 Hz). From 3e-5 it adds
  regular **doublets** (2 spikes, 20-30 ms apart, every 140-200 ms), which give CV 0.7-0.9. This is
  regular patterned firing, not burst/pause.
- **The Poisson-surprise rule cannot call short, regular groups bursts.** Even the reproduced
  triplet bursting of the old cell has burst share 0: three spikes in ~26 ms at a 16 Hz mean rate have
  a surprise of ~2, below S = 3. The rule was fixed beforehand and is not changed here. With the
  >= 3-spike burst definition of `single_cell.py`, no Kir 0.3 cell bursts regularly either. Both
  definitions therefore agree on the Kir 0.3 map.
- **Network screen: not run.** No admissible cell exists at Kv3.1 >= 1.15, or anywhere else.
