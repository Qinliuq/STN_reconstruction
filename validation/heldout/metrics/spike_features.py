"""
Per-unit spike-train features, shared by the rat data and the model (one function, two inputs).

Port of rat_6ohda/extract_features.py (the reference implementation that produced unit_cluster_features.csv);
the definitions are kept exactly, including:
- rate = n_spikes / (last - first spike), not the recording length
- ISIs <= 0 are dropped
- pause_time_frac = summed duration of ISIs > 4x mean ISI / spike-train span
- short_isi_rate_relative = fraction of ISIs < 0.5x mean ISI
- burst_share: Poisson-surprise bursts (Legendy & Salcman 1985): cores = runs of ISIs < 0.5x mean ISI, extended one
  spike at a time left / right while the surprise S increases, kept if S >= 3; fraction of spikes in kept bursts
- log_isi_bimodality: Sarle's bimodality coefficient of log10(ISI)
- units need >= 100 spikes and >= 99 positive ISIs (otherwise unit_features returns None)
Spike times are in seconds.
"""
import numpy as np
from scipy import stats as sstats

FEATURES = ['rate_hz', 'cv', 'lv', 'frac_isi_lt10ms', 'short_isi_rate_relative',
            'pause_time_frac', 'burst_share', 'log_isi_bimodality']
MIN_SPIKES = 100
MIN_ISIS = 99


def _surprise(n_spk, T, rate):
    if T <= 0 or n_spk <= 0:
        return 0.0
    p = sstats.poisson.sf(n_spk - 1, rate * T)
    return -np.log10(max(p, 1e-300))


def poisson_surprise_bursts(spike_times, isi, rate, s_min=3.0):
    """Boolean mask of spikes inside kept Poisson-surprise bursts (reference algorithm, see module docstring)."""
    n = len(spike_times)
    mask = np.zeros(n, dtype=bool)
    if n < 5 or rate <= 0:
        return mask
    in_core = isi < (1.0 / rate) / 2.0
    cores, i, m = [], 0, len(in_core)
    while i < m:
        if in_core[i]:
            j = i
            while j < m and in_core[j]:
                j += 1
            cores.append((i, j))
            i = j
        else:
            i += 1
    for s, e in cores:
        S = _surprise(e - s + 1, spike_times[e] - spike_times[s], rate)
        lo, hi = s, e
        improved = True
        while improved:
            improved = False
            if lo > 0:
                S2 = _surprise(hi - (lo - 1) + 1, spike_times[hi] - spike_times[lo - 1], rate)
                if S2 > S:
                    S, lo, improved = S2, lo - 1, True
            if hi < n - 1:
                S2 = _surprise((hi + 1) - lo + 1, spike_times[hi + 1] - spike_times[lo], rate)
                if S2 > S:
                    S, hi, improved = S2, hi + 1, True
        if S >= s_min:
            mask[lo:hi + 1] = True
    return mask


def bimodality_coeff(x):
    x = np.asarray(x)
    n = len(x)
    if n < 5:
        return np.nan
    g = sstats.skew(x)
    k = sstats.kurtosis(x, fisher=False)
    return (g ** 2 + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))


def unit_features(spike_times):
    """The 8 features of one spike train (s), or None if it fails the spike-count criteria."""
    arr = np.sort(np.asarray(spike_times, dtype=np.float64).ravel())
    if arr.size < MIN_SPIKES:
        return None
    isi = np.diff(arr)
    isi = isi[isi > 0]
    if isi.size < MIN_ISIS:
        return None
    T = arr[-1] - arr[0]
    rate = arr.size / T if T > 0 else np.nan
    mean_isi = isi.mean()
    return dict(
        n_spikes=int(arr.size), duration_s=T, rate_hz=rate,
        cv=isi.std() / mean_isi,
        lv=np.mean(3 * (isi[:-1] - isi[1:]) ** 2 / (isi[:-1] + isi[1:]) ** 2),
        frac_isi_lt10ms=np.mean(isi < 0.010),
        short_isi_rate_relative=np.mean(isi < 0.5 * mean_isi),
        pause_time_frac=isi[isi > 4 * mean_isi].sum() / T if T > 0 else np.nan,
        burst_share=poisson_surprise_bursts(arr, isi, rate).mean(),
        log_isi_bimodality=bimodality_coeff(np.log10(isi)),
    )


def features_table(trains):
    """trains: dict {unit id: spike times (s)} -> pandas DataFrame (one row per unit that passes the criteria)."""
    import pandas as pd
    rows = []
    for uid, st in trains.items():
        f = unit_features(st)
        if f is not None:
            rows.append(dict(unit=uid, **f))
    return pd.DataFrame(rows)
