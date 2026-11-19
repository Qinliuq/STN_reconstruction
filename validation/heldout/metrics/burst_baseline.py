"""
Generic burst/pause test for one spike train, against a Poisson baseline matched to its rate and length.

Used by the K_ATP screens (results/katp_kv31_map/criteria.md, results/katp_network/criteria.md). It uses the shared
Poisson-surprise burst detector (spike_features.poisson_surprise_bursts) directly, without the rat inclusion rule of
>= 100 spikes; the caller sets its own minimum (MIN_SPIKES here). Spike times are in seconds.
"""
import numpy as np

from .spike_features import poisson_surprise_bursts

MIN_SPIKES = 30
N_SURROGATES = 200


def burst_share(spike_times):
    """(burst share, CV, rate) of a train, with rate = n / (last - first spike) as in the reference features."""
    st = np.sort(np.asarray(spike_times, dtype=float))
    isi = np.diff(st)
    isi = isi[isi > 0]
    if len(st) < 3 or len(isi) < 2:
        return np.nan, np.nan, np.nan
    rate = len(st) / (st[-1] - st[0])
    return float(poisson_surprise_bursts(st, isi, rate).mean()), float(isi.std() / isi.mean()), float(rate)


def poisson_p95(rate, T, n=N_SURROGATES, seed=0):
    """95th percentile of the burst share of n homogeneous Poisson trains with this rate over T seconds."""
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        k = rng.poisson(rate * T)
        st = np.sort(rng.uniform(0, T, k))
        b, _, _ = burst_share(st)
        vals.append(0.0 if np.isnan(b) else b)
    return float(np.percentile(vals, 95))


def burst_pause(spike_times, T, seed=0):
    """Cell-level rule: >= MIN_SPIKES spikes, burst share > rate/length-matched Poisson p95 and CV > 1.
    Returns a dict with the numbers and the verdict."""
    st = np.asarray(spike_times, dtype=float)
    out = dict(n_spikes=int(len(st)), rate_hz=len(st) / T)
    if len(st) < MIN_SPIKES:
        return dict(out, burst_share=np.nan, cv=np.nan, p95_poisson=np.nan, burst_pause=False)
    b, cv, _ = burst_share(st)
    p95 = poisson_p95(len(st) / T, T, seed=seed)
    return dict(out, burst_share=b, cv=cv, p95_poisson=p95, burst_pause=bool(b > p95 and cv > 1))
