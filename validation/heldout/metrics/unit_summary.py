"""
H5 per-train summary, the one function applied to model cells and to primate segments.

ISI features (spike_features definitions, no 100-spike rule; trains with < MIN_SPIKES spikes return None) plus the
surrogate oscillation test (oscillation.py). Spike times in seconds; [t0, t1) is the analysis window.
"""
import numpy as np

from . import spike_features as sf
from .oscillation import oscillation_test

MIN_SPIKES = 20
FEATURES = ['rate_hz', 'cv', 'lv', 'pause_time_frac', 'frac_isi_lt10ms', 'burst_share']


def h5_summary(spike_times, t0, t1, seed=0, n_surr=200):
    st = np.sort(np.asarray(spike_times, dtype=float))
    st = st[(st >= t0) & (st < t1)]
    f = sf.unit_features(st, min_spikes=MIN_SPIKES, min_isis=MIN_SPIKES - 1)
    if f is None:
        return None
    o = oscillation_test(st, t0, t1, n_surr=n_surr, seed=seed)
    out = {k: float(f[k]) for k in FEATURES}
    out.update(n_spikes=int(len(st)), oscillatory=o['oscillatory'], peak_hz=o['peak_hz'], n_sig_windows=o['n_sig_windows'])
    out.update({f'ratio_{c:g}': float(r) for c, r in zip(o['centres'], o['ratio'])})
    return out
