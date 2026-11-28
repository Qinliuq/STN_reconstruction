"""
Spike phase-locking to a continuous reference signal, shared by data (EEG) and model (cortical input rate /
GPe beta drive): band-pass (zero-phase Butterworth), Hilbert phase, spike phases, PPC and preferred phase.

PPC = pairwise phase consistency (Vinck et al. 2010, PPC0): (|sum exp(i theta)|^2 - N) / (N (N - 1)), the mean
cosine of the phase difference over all spike pairs; unbiased by the spike count (unlike the PLV).
"""
import numpy as np
from scipy.signal import butter, hilbert, sosfiltfilt, welch

BANDS = {'slow': (0.5, 2.0), 'beta': (13.0, 30.0)}


def bandpass(x, fs, lo, hi, order=2):
    sos = butter(order, [lo, hi], btype='bandpass', fs=fs, output='sos')
    return sosfiltfilt(sos, np.asarray(x, dtype=float))


def band_phase(x, fs, band):
    lo, hi = BANDS[band] if isinstance(band, str) else band
    return np.angle(hilbert(bandpass(x, fs, lo, hi)))


def spike_phases(spike_times, phase, fs, t0=0.0, mask=None):
    """Phase (rad) of `phase` (sampled at fs from t0) at each spike; spikes outside the signal, or where the
    boolean sample mask is False, are dropped."""
    idx = np.round((np.asarray(spike_times) - t0) * fs).astype(int)
    ok = (idx >= 0) & (idx < len(phase))
    idx = idx[ok]
    if mask is not None:
        idx = idx[mask[idx]]
    return phase[idx]


def ppc(theta):
    theta = np.asarray(theta)
    n = len(theta)
    if n < 2:
        return np.nan
    s = np.abs(np.exp(1j * theta).sum()) ** 2
    return float((s - n) / (n * (n - 1)))


def preferred_phase(theta):
    return float(np.angle(np.exp(1j * np.asarray(theta)).mean())) if len(theta) else np.nan


def rayleigh_p(theta):
    n = len(theta)
    if n < 2:
        return np.nan
    R = np.abs(np.exp(1j * np.asarray(theta)).sum())
    return float(np.exp(np.sqrt(1 + 4 * n + 4 * (n ** 2 - R ** 2)) - (1 + 2 * n))) # Zar 1999 approximation


def locking(spike_times, ref, fs, band, t0=0.0, mask=None, min_spikes=50):
    """PPC, preferred phase, Rayleigh p and spike count of one unit to one reference signal in one band."""
    th = spike_phases(spike_times, band_phase(ref, fs, band), fs, t0, mask)
    if len(th) < min_spikes:
        return dict(n=len(th), ppc=np.nan, pref_phase=np.nan, rayleigh_p=np.nan)
    return dict(n=len(th), ppc=ppc(th), pref_phase=preferred_phase(th), rayleigh_p=rayleigh_p(th))


def slow_fraction_windows(x, fs, win_s=4.0, slow=(0.5, 2.0), total=(0.5, 45.0)):
    """Slow-wave power fraction (slow / total band) in consecutive windows; returns (window starts in s, fraction)."""
    n = int(win_s * fs)
    starts, frac = [], []
    for i in range(0, len(x) - n + 1, n):
        f, P = welch(x[i:i + n], fs=fs, nperseg=min(n, int(2 * fs)))
        band = lambda lo, hi: P[(f >= lo) & (f <= hi)].sum()
        starts.append(i / fs)
        frac.append(band(*slow) / band(*total) if band(*total) > 0 else np.nan)
    return np.array(starts), np.array(frac)


def window_mask(n_samples, fs, starts, keep, win_s=4.0):
    m = np.zeros(n_samples, dtype=bool)
    for s, k in zip(starts, keep):
        if k:
            m[int(s * fs):int((s + win_s) * fs)] = True
    return m
