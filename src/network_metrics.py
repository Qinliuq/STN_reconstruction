"""
Summary metrics of STN network activity, written after each simulation (init.py) to
<cfg.filename>_metrics.json (and the STN population-rate PSD to <cfg.filename>_psd.npz).

All metrics use t >= cfg.metrics_t0 (ms) to skip the initial transient:
- rate_PVP_pop / rate_PVN_pop: mean firing rate (Hz)
- frac_spikes_in_bursts: fraction of STN spikes in bursts (>= 3 spikes with ISI < 20 ms,
  flanked by pauses > 3x the intra-burst ISI; same detector as single_cell.py)
- beta_pow / low_pow_3_12: power of the STN population rate (1 ms bins, Welch) in 13-30 Hz / 3-12 Hz
- beta_frac / low_frac_3_12: the same, as a fraction of 1-100 Hz power
- peak_hz: frequency of the PSD maximum in 1-50 Hz
- sync_corr10ms: mean pairwise correlation of STN spike counts in 10 ms bins
- gpe_rate: mean GPe group rate (Hz)
"""
import json

import numpy as np
from scipy.signal import welch

STN_POPS = ['PVP_pop', 'PVN_pop']


def find_bursts(spikes, max_isi=20.0, min_spikes=3, pause_ratio=3.0):
    isi = np.diff(spikes)
    runs, run = [], [0]
    for i, x in enumerate(isi):
        if x < max_isi:
            run.append(i + 1)
        else:
            runs.append(run)
            run = [i + 1]
    runs.append(run)
    bursts = []
    for run in runs:
        if len(run) < min_spikes or run[0] == 0 or run[-1] == len(spikes) - 1:
            continue
        intra = np.mean(np.diff(spikes[run]))
        if min(isi[run[0] - 1], isi[run[-1]]) > pause_ratio * intra:
            bursts.append(spikes[run])
    return bursts


def compute_metrics(sim, cfg, gpe=None):
    t0, t1 = cfg.metrics_t0, cfg.duration
    T = (t1 - t0) / 1000
    spk = np.array(sim.allSimData['spkid'])
    spkt = np.array(sim.allSimData['spkt'])
    gidsByPop = {pop: list(sim.net.pops[pop].cellGids) for pop in STN_POPS}
    stnGids = [g for pop in STN_POPS for g in gidsByPop[pop]]
    m = {}
    for pop, gids in gidsByPop.items():
        m[f'rate_{pop}'] = float((np.isin(spk, gids) & (spkt >= t0)).sum() / len(gids) / T)

    inBursts = total = 0
    for g in stnGids:
        s = np.sort(spkt[spk == g])
        inBursts += sum(len(b) for b in find_bursts(s) if b[0] >= t0)
        total += (s >= t0).sum()
    m['frac_spikes_in_bursts'] = float(inBursts / total) if total else 0.0

    sel = np.isin(spk, stnGids) & (spkt >= t0)
    hist = np.histogram(spkt[sel], bins=np.arange(t0, t1 + 1, 1.0))[0].astype(float)
    f, P = welch(hist - hist.mean(), fs=1000, nperseg=min(512, len(hist)))
    band = lambda lo, hi: P[(f >= lo) & (f <= hi)].sum()
    total = band(1, 100)
    m['beta_pow'], m['low_pow_3_12'] = float(band(13, 30)), float(band(3, 12))
    m['beta_frac'], m['low_frac_3_12'] = float(band(13, 30) / total), float(band(3, 12) / total)
    inRange = (f >= 1) & (f <= 50)
    m['peak_hz'] = float(f[inRange][np.argmax(P[inRange])])

    counts = np.array([np.histogram(spkt[spk == g], bins=np.arange(t0, t1 + 1, 10.0))[0] for g in stnGids], float)
    counts = counts[counts.std(1) > 0]
    if len(counts) > 1:
        R = np.corrcoef(counts)
        m['sync_corr10ms'] = float(R[np.triu_indices(len(R), 1)].mean())
    else:
        m['sync_corr10ms'] = float('nan')

    if gpe is not None:
        rates = np.array(gpe.rates)
        m['gpe_rate'] = float(rates[int(t0 / gpe.dt):].mean())
    return m, f, P


def save_metrics(sim, cfg, gpe=None):
    m, f, P = compute_metrics(sim, cfg, gpe)
    with open(f'{cfg.filename}_metrics.json', 'w') as fh:
        json.dump(m, fh, indent=1)
    np.savez(f'{cfg.filename}_psd.npz', f=f, P=P)
    print('Metrics:', ', '.join(f'{k} {v:.4g}' for k, v in m.items()))
    return m
