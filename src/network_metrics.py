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
- profile: the same quantities along the STN axis (axis = (xnorm + ynorm) / 2, 0..1), computed at
  positions 0.05..0.95 with Gaussian weights (sigma = PROFILE_KERNEL) over the cells' positions: beta power
  shared between cells (weighted mean pairwise beta cross-power, independent of local cell count;
  beta_pow_raw = beta power of the weighted mean spike train), beta fraction of the weighted mean train, weighted mean pairwise synchrony, weighted mean per-cell
  burst fraction and rate, and the effective number of cells per position (n_eff)
- beta_peak_axis / sync_peak_axis: positions of the maxima of the beta-power and synchrony profiles
  (among positions with n_eff >= PROFILE_MIN_CELLS)
"""
import json

import numpy as np
from scipy.signal import welch

STN_POPS = ['PVP_pop', 'PVN_pop']
PROFILE_POSITIONS = np.round(np.arange(0.05, 0.951, 0.05), 2)
PROFILE_KERNEL = 0.12
PROFILE_MIN_CELLS = 5 # minimum n_eff for a position to count for the peak


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
    cellBurstFrac, cellRate = [], []
    for g in stnGids:
        s = np.sort(spkt[spk == g])
        nb = sum(len(b) for b in find_bursts(s) if b[0] >= t0)
        ns = (s >= t0).sum()
        inBursts += nb
        total += ns
        cellBurstFrac.append(nb / ns if ns else 0.0)
        cellRate.append(ns / T)
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
    active = counts.std(1) > 0
    if active.sum() > 1:
        R = np.corrcoef(counts[active])
        m['sync_corr10ms'] = float(R[np.triu_indices(len(R), 1)].mean())
    else:
        R = None
        m['sync_corr10ms'] = float('nan')

    # profiles along the STN axis
    cellsByGid = {cell.gid: cell for cell in sim.net.cells}
    axis = np.array([np.clip((cellsByGid[g].tags['xnorm'] + cellsByGid[g].tags['ynorm']) / 2, 0, 1) for g in stnGids])
    trains = np.array([np.histogram(spkt[spk == g], bins=np.arange(t0, t1 + 1, 1.0))[0] for g in stnGids], float)
    nperseg = min(512, trains.shape[1])
    fc, Pc = welch(trains - trains.mean(1, keepdims=True), fs=1000, nperseg=nperseg, axis=1)
    cellBeta = Pc[:, (fc >= 13) & (fc <= 30)].sum(1)  # each cell's own beta power
    prof = {k: [] for k in ('beta_pow', 'beta_pow_raw', 'beta_frac', 'sync', 'bursts', 'rate', 'n_eff')}
    for x in PROFILE_POSITIONS:
        w = np.exp(-(axis - x) ** 2 / (2 * PROFILE_KERNEL ** 2))
        prof['n_eff'].append(float(w.sum()))
        wn = w / w.sum()
        trace = wn @ trains
        fx, Px = welch(trace - trace.mean(), fs=1000, nperseg=nperseg)
        bandx = lambda lo, hi: Px[(fx >= lo) & (fx <= hi)].sum()
        # beta power shared between cells: power of the weighted mean minus the cells' own (auto) terms, normalized
        # by the pair weights = weighted mean pairwise beta cross-power, independent of how many cells are nearby
        prof['beta_pow'].append(float((bandx(13, 30) - (wn ** 2) @ cellBeta) / (1 - (wn ** 2).sum())))
        prof['beta_pow_raw'].append(float(bandx(13, 30)))
        prof['beta_frac'].append(float(bandx(13, 30) / bandx(1, 100)) if bandx(1, 100) > 0 else 0.0)
        prof['bursts'].append(float(wn @ np.array(cellBurstFrac)))
        prof['rate'].append(float(wn @ np.array(cellRate)))
        if R is not None:
            wa = w[active]
            W = np.outer(wa, wa)
            iu = np.triu_indices(len(wa), 1)
            prof['sync'].append(float((W[iu] * R[iu]).sum() / W[iu].sum()))
        else:
            prof['sync'].append(float('nan'))
    m['profile'] = {'positions': PROFILE_POSITIONS.tolist(), **prof}
    # peaks only where the window covers enough cells: the weighted mean spike train of few cells is noisier,
    # which inflates its power at the sparsely populated ends of the axis
    ok = np.array(prof['n_eff']) >= PROFILE_MIN_CELLS
    peak = lambda k: float(PROFILE_POSITIONS[ok][int(np.nanargmax(np.array(prof[k])[ok]))]) if ok.any() else float('nan')
    m['beta_peak_axis'] = peak('beta_pow')
    m['sync_peak_axis'] = peak('sync') if R is not None else float('nan')

    if gpe is not None:
        rates = np.array(gpe.rates)
        m['gpe_rate'] = float(rates[int(t0 / gpe.dt):].mean())
    return m, f, P


def save_metrics(sim, cfg, gpe=None):
    m, f, P = compute_metrics(sim, cfg, gpe)
    with open(f'{cfg.filename}_metrics.json', 'w') as fh:
        json.dump(m, fh, indent=1)
    np.savez(f'{cfg.filename}_psd.npz', f=f, P=P)
    print('Metrics:', ', '.join(f'{k} {v:.4g}' for k, v in m.items() if not isinstance(v, dict)))
    return m
