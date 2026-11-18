# Kv3.1 x K_ATP single-cell map (PV+): admissibility criteria (fixed before the runs)

Written 2026-11-18, before any run of this map. The cells are isolated PV+ cells, k_Ca = 20 uM.
The protocols are those of `src/single_cell.py`, run for 11 s (`--duration 11000`; statistics over
1-11 s), plus the 4 s depolarizing set (`--protocols depol`). There are no rat values here.

A (Kv3.1, K_ATP gbar) cell is **admissible** when all three criteria hold.

1. **Spontaneous burst/pause** (cell-level rule of `results/katp_network/criteria.md`):
   - at least 30 spikes over 1-11 s;
   - Poisson-surprise burst share above the 95th percentile of 200 homogeneous Poisson trains
     with the same rate and length;
   - ISI CV above 1.
2. **No depolarization block at +0.8 nA, 4 s** (protocol of `results/kv31_single_cell/`):
   - the cell still spikes in the last 500 ms;
   - the trough over the last second is below -50 mV, i.e. fully repolarized spikes rather than
     the near-block oscillation (trough ~ -35 mV).
3. **Phenotype preserved**:
   - at -0.16 nA: rhythmic bursting, i.e. at least 1 burst/s over 1-11 s with a mean of at least
     3 spikes per burst (burst = >= 3 spikes with ISI < 20 ms, flanked by pauses > 3x the
     intra-burst ISI; `single_cell.find_bursts`);
   - rebound after -0.25 nA (1000-1500 ms): at least 8 spikes, at a mean frequency of at least
     100 Hz.

The diagnostic cells (Kv3.1 1.0 / Kir 0.5, Kv3.1 1.0 / Kir 0.3, Kv3.1 1.2 / Kir 0.5) are reported
with the same measures. The reproduction check is against the 2026-10-18 result: ~5 Hz bursts of
~3 spikes, CV ~1.1, at 3e-5 with Kv3.1 1.0 / Kir 0.5. It uses the 1-4 s window of the original run,
reported alongside the 1-11 s numbers.
