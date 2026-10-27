# K_ATP (katp_stn.mod) in isolated PV+ and PV- cells, 2026-10-18

Ca2+-gated K_ATP (Kir6.2/SUR1; `mod/katp_stn.mod`), inserted uniformly in soma and dendrites
(`tools.add_katp`), swept with `src/single_cell.py --param katp_gbar` (model state: commit `ac8cd07`;
other conductances as in `src/cfg.py`, pdscale = 1, ATP/ADP defaults so fatp ~ 0.53). Each line in
`*.txt` gives rate, ISI CV and bursting over 1-4 s, plus the mean somatic [Ca]i and K_ATP open
probability `o`. Figures: voltage traces per protocol (rows) and gbar (columns).

```
python src/single_cell.py --cell PVP --param katp_gbar --values 0 3e-5 1e-4 3e-4 1e-3
python src/single_cell.py --cell PVP --param katp_gbar --values 0 1e-5 3e-5 1e-4 3e-4 --cfg katp_kca=0.02
```

## The Ca2+ gate saturates with the mod default kca

This model's Ca pool (`Cacum`, 200 nm shell) reaches ~20-30 uM during spontaneous firing and
70-80 uM at +0.05 nA, and stays at 0.2-1.2 uM in silent cells. With the mod default
kca = 0.5 uM, K_ATP is ~fully open (o ~ 0.99) in any firing cell and half open even in silent ones,
so it acts as a near-constant K+ leak. Rescaling kca to 20 uM (`cfg.katp_kca = 0.02`, between the
silent and firing levels) restores activity dependence (o ~ 0.6 while firing spontaneously, ~0.9
at +0.05 nA, ~0.01-0.1 when silent). The default is left unchanged in `cfg.py`; kca should be
chosen to match this model's Ca scale, not the physiological 0.5 uM.

## Effect on firing (kca = 20 uM; the 0.5 uM runs are qualitatively the same, slightly stronger)

| gbar (S/cm2) | PV+ spontaneous | PV+ rebound | PV- spontaneous | PV- rebound |
|---|---|---|---|---|
| 0 | 17 Hz tonic (CV 0.01) | 14 spikes / 79 ms | 11 Hz tonic | 12 spikes / 71 ms |
| 1e-5 | 14 Hz, irregular (CV 0.6) | 14 / 82 ms | **silent** | 10 / 57 ms |
| 3e-5 | **bursting**: ~5 Hz bursts of 3 spikes (CV 1.1) | 13 / 73 ms | silent | 9 / 49 ms |
| 1e-4 | silent | 10 / 54 ms | silent | 8 / 46 ms |
| 3e-4 | silent | 7 / 34 ms | silent | 5 / 21 ms |

With depolarizing drive (+0.05 nA) both types keep firing up to 3e-5 (PV+ 55 -> 16 Hz,
PV- 50 -> 13 Hz) and are silenced from 1e-4.

- K_ATP is a slow (tau = 300 ms), activity-dependent K+ current: firing raises [Ca]i, which opens
  K_ATP and terminates firing; [Ca]i then decays and firing resumes. At low density this turns PV+
  tonic pacemaking into rhythmic bursting.
- PV- cells, weaker pacemakers, are silenced by the smallest density tested, while PV+ cells burst;
  K_ATP therefore sharpens the PV+ (bursting) / PV- (non-bursting) difference.
- At moderate density (1e-4) cells are silent at rest but still fire a rebound burst after
  hyperpolarization (PV+ 10 spikes), i.e. they fire only on release from inhibition.
- Higher densities shorten and then abolish rebound bursts and silence both types.
- Note that the density is applied to all sections (dendrites carry most of the membrane area), so
  even 1e-5 S/cm2 is a substantial total conductance.

Implications for the network: augmented K_ATP after dopamine depletion (`katp_pdscale` > 1;
Shen & Johnson 2012) is a candidate mechanism for more bursty parkinsonian STN firing, and K_ATP
opening (suppressing firing) acts in the opposite direction to Kir block.
