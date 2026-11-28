# K_ATP network screen: selection criteria (fixed before the screen ran)

Written 2026-11-18, before any network run with K_ATP. The thresholds use only model-internal
quantities, generic burst definitions and the configuration-D readouts of the results draft.
No rat feature values are used: no CV, pause fraction, burst share or cluster proportion from the
H1/H2 data. The analysis is in `results/katp_network/analyze.py`, and every time is t >= 0.5 s.

A setting passes when all four criteria hold. The pick is the **smallest PV+ gbar** that passes,
using the smaller PV- gbar on a tie. If no setting passes, that is the result; the search is not
widened.

## 1. PV+ burst/pause, PV- comparatively tonic

Per cell, from `validation/heldout/metrics/spike_features.py`: the ISI CV, and the Poisson-surprise
burst share (Legendy & Salcman, S >= 3). The screen applies the functions directly, without the
100-spike inclusion rule, but requires >= 30 spikes per cell.

The Poisson baseline is matched to each cell's rate and train length: 200 homogeneous Poisson
trains, and their 95th percentile of burst share (`p95_poisson`).

A cell is **burst/pause** when its burst share is above its own `p95_poisson` **and** its CV is
above 1.

Pooled over seeds, the setting must have:
- (1a) a burst/pause fraction of at least 0.5 among PV+ cells;
- (1b) a PV+ median CV above 1, and a PV+ median burst share above the PV+ median `p95_poisson`;
- (1c) PV- comparatively tonic: PV- median CV below the PV+ median CV, and a PV- burst/pause
  fraction below half the PV+ fraction.

## 2. Plausible PD rates

The seed-mean population rates must be 15-45 Hz for PV+ and 8-30 Hz for PV-. No cell may be
silent (< 1 Hz over the analysed window) in any seed.

## 3. PD beta features of configuration D still hold

- (3a) The beta fraction of the population rate is above the healthy configuration-D value. Seed
  means are compared, and the healthy runs use the same seeds and duration.
- (3b) The peak of the seed-mean beta-power axis profile (shared beta) lies at axis 0.55-0.75.
- (3c) That profile at 0.65 is larger than at 0.25 (the upper-middle maximum, not a flat profile).

## 4. Healthy state unchanged

The healthy K_ATP density is 0 (`cfg.katp_state_scale = {'healthy': 0, 'parkinsonian': 1}`), so the
channel is not inserted in healthy runs. This is verified once, by bit-identical spikes for a
healthy run with a nonzero `katp_gbar` against one without.
