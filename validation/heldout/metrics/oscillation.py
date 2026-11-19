"""
Single-train oscillation test against ISI-shuffled surrogates (one function for model and data trains).

Written for the dropped H5 test (see validation/heldout/README.md); kept as a tested utility. Parameters:
- 1 ms bins, mean-subtracted; Welch PSD (Hann, nperseg 2048 bins, 50% overlap)
- 2 Hz windows with edges 3, 5, ..., 41 Hz; window power = mean PSD in [lo, hi)
- 200 surrogates: the same ISIs in random order, starting at the first spike
- oscillatory if any window's power exceeds the 99.9th percentile of that window's surrogate powers
- peak = centre of the significant window with the largest power / surrogate-median ratio
Spike times in seconds.
"""
import numpy as np
from scipy.signal import welch

BIN_S = 0.001
NPERSEG = 2048
EDGES = np.arange(3, 42, 2)
N_SURR = 200
PCT = 99.9


def window_power(spike_times, t0, t1):
    n = int(round((t1 - t0) / BIN_S))
    x = np.zeros(n)
    idx = ((np.asarray(spike_times) - t0) / BIN_S).astype(int)
    idx = idx[(idx >= 0) & (idx < n)]
    np.add.at(x, idx, 1.0)
    f, P = welch(x - x.mean(), fs=1 / BIN_S, window='hann', nperseg=min(NPERSEG, n), noverlap=min(NPERSEG, n) // 2)
    return np.array([P[(f >= lo) & (f < hi)].mean() for lo, hi in zip(EDGES[:-1], EDGES[1:])])


def oscillation_test(spike_times, t0, t1, n_surr=N_SURR, seed=0):
    """dict: oscillatory (bool), peak_hz (nan if none), n_sig_windows, ratio (power / surrogate median per window),
    power, thresholds. Trains with < 3 spikes are not oscillatory."""
    st = np.sort(np.asarray(spike_times, dtype=float))
    st = st[(st >= t0) & (st < t1)]
    centres = (EDGES[:-1] + EDGES[1:]) / 2
    if len(st) < 3:
        return dict(oscillatory=False, peak_hz=np.nan, n_sig_windows=0, centres=centres,
                    ratio=np.full(len(centres), np.nan), power=np.full(len(centres), np.nan))
    p = window_power(st, t0, t1)
    rng = np.random.default_rng(seed)
    isi = np.diff(st)
    S = np.empty((n_surr, len(p)))
    for k in range(n_surr):
        S[k] = window_power(st[0] + np.concatenate([[0.0], np.cumsum(rng.permutation(isi))]), t0, t1)
    thr = np.percentile(S, PCT, axis=0)
    med = np.median(S, axis=0)
    sig = p > thr
    ratio = p / med
    peak = float(centres[sig][np.argmax(ratio[sig])]) if sig.any() else np.nan
    return dict(oscillatory=bool(sig.any()), peak_hz=peak, n_sig_windows=int(sig.sum()), centres=centres,
                ratio=ratio, power=p, thresholds=thr)
