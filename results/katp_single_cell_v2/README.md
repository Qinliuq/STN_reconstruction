# K_ATP in isolated PV+ and PV- cells, configuration D, 2026-11-18

These runs repeat `results/katp_single_cell/` (2026-10-18) on the current cell model, configuration D:
- equal Kir 0.3 / 0.3;
- Kv3.1 PV+ 1.2 / PV- 1.15.

The earlier sweep predates both changes. Here k_Ca = 20 uM (`katp_kca = 0.02` mM), pdscale = 1, and
the ATP/ADP defaults give fatp ~ 0.53. The model state is commit `95a580c`.

```
python src/single_cell.py --cell PVP --param katp_gbar --values 0 5e-6 1e-5 2e-5 3e-5 5e-5 --cfg katp_kca=0.02
python src/single_cell.py --cell PVN --param katp_gbar --values 0 1e-6 3e-6 1e-5 --cfg katp_kca=0.02
```

Files:
- `pvp_kca20uM.txt`, `pvn_kca20uM.txt`: full output. Rates and CV are over 1-4 s. A burst is >= 3
  spikes with ISI < 20 ms, flanked by pauses.
- `PVP_kca20uM.png`, `PVN_kca20uM.png`: voltage traces, one row per protocol and one column per gbar.

## Results

| gbar (S/cm2) | PV+ spontaneous | PV+ +0.05 nA | PV+ hold -0.16 nA | PV+ rebound |
|---|---|---|---|---|
| 0 | 22.3 Hz tonic (CV 0.00) | 59 Hz | 13 Hz, 3-spike bursts (CV 1.36) | 12 spikes / 69 ms |
| 5e-6 | 15.3 Hz tonic | 54 Hz | 12.7 Hz bursts | 12 / 70 ms |
| 1e-5 | 10.7 Hz tonic | 48 Hz | 11.0 Hz bursts | 12 / 71 ms |
| 2e-5 | 7.0 Hz tonic | 35 Hz | 8.7 Hz bursts | 11 / 62 ms |
| 3e-5 | 5.7 Hz tonic (CV 0.01) | 21 Hz | **silent** | 11 / 62 ms |
| 5e-5 | 5.7 Hz, skipped spikes (CV 0.33) | 7 Hz | silent | 10 / 57 ms |

| gbar (S/cm2) | PV- spontaneous | PV- +0.05 nA | PV- rebound |
|---|---|---|---|
| 0 | 6.7 Hz tonic | 40 Hz | 8 / 44 ms |
| 1e-6 | **silent** | 38 Hz | 8 / 44 ms |
| 3e-6 | silent | 36 Hz | 8 / 44 ms |
| 1e-5 | silent | 26 Hz | 8 / 45 ms |

Neither type fires at -0.25 or -0.35 nA at any density.

## Interpretation

- **K_ATP no longer makes PV+ cells burst.** In the pre-configuration-D cells, 3e-5 S/cm2 turned PV+
  pacemaking into ~5 Hz rhythmic 3-spike bursts (CV 1.1). In configuration-D cells, K_ATP only
  slows tonic pacemaking: 22 -> 6 Hz with CV <= 0.01 up to 3e-5. At 5e-5 a subthreshold oscillation
  appears and some spikes are skipped (CV 0.33), but there are still no bursts.
- **K_ATP removes the bursting that is already there.** The only PV+ bursting is under the
  -0.16 nA hold, and it is present without K_ATP. K_ATP reduces it and abolishes it from 3e-5.
- **Rebound bursts survive.** PV+ keeps a 10-12 spike rebound burst at every density.
- **PV- cells are silenced at rest by the smallest density, 1e-6.** As before, PV- are weaker
  pacemakers at 6.7 Hz. With depolarizing drive they keep firing at every density tested (26-38 Hz).
- The likely cause is the raised Kv3.1 and the equal Kir. These configuration-D changes made
  pacemaking more robust, and that robustness now dominates the slow Ca2+-gated K+ feedback.
  This explanation was not tested further, and the search was not widened (kca, tau).

Consequence for the network screen: the step-2 grid was defined as PV+ densities where single cells
burst but are not silent. In configuration D that range is empty.
