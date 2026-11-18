Accepted by the user on 2026-11-18, before any H5 simulation or primate download.

What was searched:
- Terms: all files, all git branches and the plan documents (`doc/Plan.docx`, `doc/Plan.pdf`) for Galvan, Wichmann,
  MPTP, primate, monkey, macaque, rhesus, zenodo and 6581800. There were no hits.
- Numbers: 22.9, 36.8, 1.46 and 1.32. The only hits are coincidental: rat CV 1.32 in H2, a phase of −1.32 rad in H3,
  a rescue gap of 1.32.
- `primate-isi-data-notes.md` is not in the repository or on the analysis machine.
- The transcript of the calibration session (`session_01V1zyES4ZiM8ZqYqCwKTb3X`, 2026-10-02 to 10-27) was sampled
  only (the last ~2.5 h, no hits), not read in full.

Calibration steps after 2026-10-09:

| date | commit | step | tuned to |
|---|---|---|---|
| 10-10/11 | `b88b00b`, `2471c3f` | cortical/GPe input generators (20 Hz beta, GPe 25 Hz, cortex 80 Hz) | none; labelled "illustrative defaults" |
| 10-17 | `e0d9960` | STN–STN weight, cortex 1 nS | PD distinct from healthy in beta and synchrony; inputs rate-matched, so equal rates |
| 10-17 | `246c7c0`, `a7a9946`, `8396160` | GPe as beta source, 8 ms GABA, 2 bursts/s | STN beta and synchrony ratio PD/healthy, peak ~20 Hz; reliable burst coverage |
| 10-18 | `94669cf`, `7f0fb90` | topographic GPe, cortex 0.7 nS | PD bursting ≥ 1.5× healthy; beta/sync peak at axis ~0.65; "plausible healthy rates ~30/20 Hz" |
| 10-24 | `347b84f`, `45f6e04` | Kv3.1 1.2/1.15, cortex 0.8 nS | single-cell depolarization-block protection; beta peak 0.6–0.7 |
| 10-25 | `7a57b8d`, `4e0eb7e` | equal Kir, PD synapse scale 2.0 | single-cell phenotype; restoring the PD state |
| 10-26 | `b3cf1c1` | GABA scale PV+ 1.13 / PV− 0.87 | GABA:glutamate ratio |

H5 targets:

| H5 target | class | evidence |
|---|---|---|
| Rate change healthy → PD (data ×1.6) | **clean** | PD and healthy inputs were rate-matched by design. The model's PD rate *falls* (PV+ 37.5 → 32.9, PV− 20.3 → 16.4 Hz) and was never tuned toward a rise. |
| Irregularity unchanged (CV, LV, pause fraction) | **clean** | The calibration target was the opposite: PD ≥ 1.5× burstier than healthy. |
| Oscillatory fraction (6/13 → 24/27) | **clean (numbers); direction built in** | The primate fractions were never used. "PD oscillates more" is the definition of the PD state, so the direction can't fail independently of the input assumption. |
| Oscillation frequency (MPTP peaks 4–17 Hz) | **clean** | 20 Hz was a default from 10-10, before any calibration, and contradicts the primate peaks. It tests the input assumption only. |
| Within-state levels (secondary) | **unclear** | Healthy rates were targeted as "plausible ~30/20 Hz" with no stated source. The model's healthy PV− rate of 20.3 Hz sits near the primate normal median of 22.9 Hz. This affects only secondary statistics, not the four change tests. |
