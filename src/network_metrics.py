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
- per cell type (<pop> = PVP_pop / PVN_pop), independent of the number of cells of each type:
  beta_pow_pair_<pop>: mean pairwise beta cross-power (13-30 Hz) between cells of that type (shared beta);
  beta_pow_cell_<pop>: mean beta power of the single cells' spike trains (own beta, incl. their spiking noise);
  beta_frac_<pop>: beta fraction (13-30 Hz / 1-100 Hz) of the type's mean spike train;
  sync_<pop>: mean pairwise correlation of 10 ms counts within the type.
  beta_pow_pair_STN / beta_pow_pair_cross: the same pairwise beta cross-power over all STN pairs / PV+-PV- pairs.

Spike times are saved with the metrics (<filename>_spikes.npz), so metrics can be recomputed after changes here:
    python src/network_metrics.py sweeps/out/<spec>/<run dir> [...]
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


def spike_data(sim, cfg, gpe=None):
    # everything the metrics need, as arrays (also what is saved to <filename>_spikes.npz)
    gidsByPop = {pop: list(sim.net.pops[pop].cellGids) for pop in STN_POPS}
    stnGids = [g for pop in STN_POPS for g in gidsByPop[pop]]
    cellsByGid = {cell.gid: cell for cell in sim.net.cells}
    spk = np.array(sim.allSimData['spkid'])
    spkt = np.array(sim.allSimData['spkt'])
    sel = np.isin(spk, stnGids)
    d = dict(spkid=spk[sel].astype(np.int32), spkt=spkt[sel], t0=cfg.metrics_t0, t1=cfg.duration,
             axis=np.array([np.clip((cellsByGid[g].tags['xnorm'] + cellsByGid[g].tags['ynorm']) / 2, 0, 1) for g in stnGids]),
             **{f'gids_{pop}': np.array(gidsByPop[pop], dtype=np.int32) for pop in STN_POPS})
    if gpe is not None:
        d['gpe_rate'] = np.array(gpe.rates).mean(1) # mean over GPe groups, one value per step of gpe_dt_ms
        d['gpe_dt_ms'] = gpe.dt
        if getattr(gpe, 'betaPhase', None):
            d['gpe_beta_phase'] = np.array(gpe.betaPhase, dtype=np.float32) # per step of gpe_dt_ms, from t = 0
            d['gpe_in_burst'] = np.array(gpe.inBurst, dtype=bool)
        at, ac = np.asarray(gpe.arrivalTimes), np.asarray(gpe.arrivalCells, dtype=int)
        d['gpe_arrivals_cell'] = np.bincount(ac[at >= cfg.metrics_t0], minlength=len(stnGids)).astype(np.int32)  # delivered (log)
    try: # cortical input signals (1 ms bins; see stn_inputs.CTX_SIGNALS)
        from stn_inputs import CTX_SIGNALS
        d.update({k: np.asarray(v) for k, v in CTX_SIGNALS.items()})
        from stn_inputs import CELL_SCALES # per-cell receptor-profile scales (only when cfg.ctx_gscale / gpe_gscale are set)
        d.update({k: np.asarray(v) for k, v in CELL_SCALES.items()})
    except ImportError:
        pass
    return d


def compute_metrics(d):
    t0, t1 = float(d['t0']), float(d['t1'])
    T = (t1 - t0) / 1000
    spk, spkt = d['spkid'], d['spkt']
    gidsByPop = {pop: list(d[f'gids_{pop}']) for pop in STN_POPS}
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
    axis = np.asarray(d['axis'])
    trains = np.array([np.histogram(spkt[spk == g], bins=np.arange(t0, t1 + 1, 1.0))[0] for g in stnGids], float)
    nperseg = min(512, trains.shape[1])
    fc, Pc = welch(trains - trains.mean(1, keepdims=True), fs=1000, nperseg=nperseg, axis=1)
    cellBeta = Pc[:, (fc >= 13) & (fc <= 30)].sum(1)  # each cell's own beta power

    # per cell type: shared (pairwise) beta, single-cell beta, beta fraction of the type's mean train, synchrony.
    # Welch's estimate is bilinear, so the power of a summed train = sum of all auto- and cross-spectra, and the mean
    # pairwise cross-power follows from the power of the sum minus the auto terms
    betaBand = lambda x: welch(x - x.mean(), fs=1000, nperseg=nperseg)
    def band_power(x, lo, hi):
        f, P = betaBand(x)
        return P[(f >= lo) & (f <= hi)].sum()
    idxByPop = {pop: np.arange(len(gidsByPop[pop])) + sum(len(gidsByPop[p]) for p in STN_POPS[:k]) for k, pop in enumerate(STN_POPS)}
    sumBeta = {}
    for pop, idx in idxByPop.items():
        n = len(idx)
        sumTrain = trains[idx].sum(0)
        sumBeta[pop] = band_power(sumTrain, 13, 30)
        m[f'beta_pow_pair_{pop}'] = float((sumBeta[pop] - cellBeta[idx].sum()) / (n * (n - 1))) if n > 1 else float('nan')
        m[f'beta_pow_cell_{pop}'] = float(cellBeta[idx].mean())
        total = band_power(sumTrain, 1, 100)
        m[f'beta_frac_{pop}'] = float(sumBeta[pop] / total) if total > 0 else 0.0
        if R is not None:
            act = np.flatnonzero(active)
            sub = np.flatnonzero(np.isin(act, idx))
            Rs = R[np.ix_(sub, sub)]
            m[f'sync_{pop}'] = float(Rs[np.triu_indices(len(sub), 1)].mean()) if len(sub) > 1 else float('nan')
        else:
            m[f'sync_{pop}'] = float('nan')
    nAll = len(stnGids)
    allBeta = band_power(trains.sum(0), 13, 30)
    m['beta_pow_pair_STN'] = float((allBeta - cellBeta.sum()) / (nAll * (nAll - 1)))
    nP, nN = (len(idxByPop[p]) for p in STN_POPS)
    m['beta_pow_pair_cross'] = float((allBeta - sumBeta['PVP_pop'] - sumBeta['PVN_pop']) / (2 * nP * nN))
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

    if 'gpe_beta_phase' in d and 'gpe_in_burst' in d:
        m['beta_ppc'] = beta_ppc(d)
    if 'gpe_rate' in d:
        m['gpe_rate'] = float(np.asarray(d['gpe_rate'])[int(t0 / float(d['gpe_dt_ms'])):].mean())
    return m, f, P


PPC_MIN_SPIKES = 10  # minimum spikes inside beta episodes for a cell's PPC
PPC_SHUFFLES = 200
PPC_MIN_SHIFT_MS = 200


def beta_ppc(d, n_shuffle=PPC_SHUFFLES, seed=0):
    """Per-cell locking to the IMPOSED pallidal beta (results/edge_artefacts, Step 3).

    d: one run's spike data, or a list of runs of the SAME network (same cells and positions, e.g. input seeds), whose
    in-episode spikes are then pooled per cell (PPC across runs). Phase = the GPe generator's common beta phase
    (d['gpe_beta_phase'], one value per GPe step), only inside beta episodes (d['gpe_in_burst']), t >= t0. Per cell:
    pairwise phase consistency over pairs of in-episode spikes from DIFFERENT beta episodes (PPC_1 of Vinck et al.
    2012): PPC = (|sum_all|^2 - sum_e |sum_e|^2) / (N^2 - sum_e n_e^2), sum = sum of e^(i phi) (all spikes / episode e).
    Unbiased by spike count and by within-episode spike-history dependence (regular, pacemaking firing within an
    episode biases the all-pairs PPC negative). NaN if N < PPC_MIN_SPIKES or fewer than 2 episodes have spikes. Significance: each run's phase and episode series are rotated together by
    a random circular shift (>= PPC_MIN_SHIFT_MS from 0), n_shuffle times; p = (1 + #(PPC_shuffled >= PPC)) /
    (1 + n_shuffle); locked if p < 0.05. Profile along the axis with the same Gaussian kernel as the other profiles,
    n_eff = sum of kernel weights over cells with a valid PPC (cells, not pairs).
    """
    runs = d if isinstance(d, (list, tuple)) else [d]
    rng = np.random.default_rng(seed)
    gids0 = [g for pop in STN_POPS for g in runs[0][f'gids_{pop}']]
    nC, nP = len(gids0), len(runs[0]['gids_PVP_pop'])
    prep = []
    for r in runs:
        t0, dt = float(r['t0']), float(r['gpe_dt_ms'])
        phase, inb = np.asarray(r['gpe_beta_phase'], float), np.asarray(r['gpe_in_burst'], bool)
        i0 = int(np.ceil(t0 / dt))
        L = len(phase) - i0
        gids = [g for pop in STN_POPS for g in r[f'gids_{pop}']]
        spk, spkt = r['spkid'], r['spkt']
        idx = [np.clip(np.floor(spkt[(spk == g) & (spkt >= t0)] / dt).astype(int) - i0, 0, L - 1) for g in gids]
        shifts = np.concatenate([[0], rng.integers(PPC_MIN_SHIFT_MS / dt, L - PPC_MIN_SHIFT_MS / dt, n_shuffle)])
        ib = inb[i0:]
        ep = np.cumsum(np.diff(np.concatenate([[0], ib.astype(int)])) == 1)  # episode number (valid where ib)
        if ib[0] and ib[-1]:  # an episode wrapping around the circular shift is one episode
            ep[ep == ep[-1]] = ep[0]
        prep.append((phase[i0:], ib, ep, L, idx, shifts))

    def stats(c, s):  # cell c, shuffle s: resultant, sum over episodes of |episode resultant|^2, spikes, sum n_e^2
        z, se, n, sn = 0j, 0.0, 0, 0
        for r, (ph, ib, ep, L, idx, shifts) in enumerate(prep):
            j = (idx[c] + shifts[s]) % L
            j = j[ib[j]]
            if not len(j):
                continue
            e = np.unique(ep[j], return_inverse=True)[1]
            u = np.exp(1j * ph[j])
            ze = np.bincount(e, u.real) + 1j * np.bincount(e, u.imag)
            ne = np.bincount(e)
            z += u.sum()
            se += float((np.abs(ze) ** 2).sum())
            n += len(j)
            sn += int((ne ** 2).sum())
        return z, se, n, sn

    def ppc_f(z, se, n, sn):
        if n < PPC_MIN_SPIKES or n * n - sn <= 0:
            return np.nan
        return float((abs(z) ** 2 - se) / (n * n - sn))

    ppc, nsp, pval = np.full(nC, np.nan), np.zeros(nC, int), np.full(nC, np.nan)
    for c in range(nC):
        st = stats(c, 0)
        ppc[c], nsp[c] = ppc_f(*st), st[2]
        if np.isfinite(ppc[c]):
            sh = np.array([ppc_f(*stats(c, s)) for s in range(1, n_shuffle + 1)])
            sh = sh[np.isfinite(sh)]
            pval[c] = (1 + (sh >= ppc[c]).sum()) / (1 + len(sh))
    valid = np.isfinite(ppc)
    locked = valid & (pval < 0.05)
    axis = np.asarray(runs[0]['axis'])
    prof, neff = [], []
    for x in PROFILE_POSITIONS:
        w = np.exp(-(axis - x) ** 2 / (2 * PROFILE_KERNEL ** 2)) * valid
        neff.append(float(w.sum()))
        prof.append(float((w * np.nan_to_num(ppc)).sum() / w.sum()) if w.sum() > 0 else float('nan'))
    isP = np.arange(nC) < nP
    frac = lambda m: float(locked[m & valid].sum() / (m & valid).sum()) if (m & valid).any() else float('nan')
    mean = lambda m: float(np.nanmean(ppc[m])) if (m & valid).any() else float('nan')
    return dict(ppc_cell=[float(v) for v in ppc], ppc_n_spikes=nsp.tolist(), ppc_p=[float(v) for v in pval],
                ppc_mean=mean(np.ones(nC, bool)), ppc_PVP=mean(isP), ppc_PVN=mean(~isP),
                frac_locked=frac(np.ones(nC, bool)), frac_locked_PVP=frac(isP), frac_locked_PVN=frac(~isP),
                n_valid=int(valid.sum()), profile_ppc=prof, profile_ppc_n_eff=neff)


def write_metrics(d, filename):
    m, f, P = compute_metrics(d)
    with open(f'{filename}_metrics.json', 'w') as fh:
        json.dump(m, fh, indent=1)
    np.savez(f'{filename}_psd.npz', f=f, P=P)
    print('Metrics:', ', '.join(f'{k} {v:.4g}' for k, v in m.items() if not isinstance(v, dict)))
    return m


def save_metrics(sim, cfg, gpe=None):
    d = spike_data(sim, cfg, gpe)
    np.savez_compressed(f'{cfg.filename}_spikes.npz', **d)
    return write_metrics(d, cfg.filename)


if __name__ == '__main__':
    # recompute the metrics of saved runs: python src/network_metrics.py <run dir or <filename>_spikes.npz> [...]
    import os
    import sys
    for arg in sys.argv[1:]:
        path = os.path.join(arg, 'model_output_spikes.npz') if os.path.isdir(arg) else arg
        print(path)
        write_metrics(dict(np.load(path)), path[:-len('_spikes.npz')])
