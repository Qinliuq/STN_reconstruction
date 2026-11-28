"""
Responses to a repeated stimulation protocol (pre / on / post blocks), shared by data and model:
per-cycle window rates, per-unit Wilcoxon (stim vs pre across cycles, p < 0.01 as in the ChR2 analysis),
responder labels, modulation index and trial-averaged rate profiles normalized to the pre-stimulus rate.
"""
import numpy as np
from scipy.stats import wilcoxon

ALPHA = 0.01


def window_rates(spike_times, onsets, t_from, t_to):
    """Rate (Hz) in [onset + t_from, onset + t_to) for each onset."""
    st = np.sort(np.asarray(spike_times))
    a = np.searchsorted(st, np.asarray(onsets) + t_from)
    b = np.searchsorted(st, np.asarray(onsets) + t_to)
    return (b - a) / (t_to - t_from)


def classify(pre, stim, alpha=ALPHA):
    """Wilcoxon signed-rank test of stim vs pre rates across cycles -> dict(p, resp in {'exc', 'inh', 'ns'}, mi, ratio)."""
    pre, stim = np.asarray(pre, float), np.asarray(stim, float)
    d = stim - pre
    p = 1.0 if np.all(d == 0) else float(wilcoxon(stim, pre).pvalue)
    resp = 'ns' if p >= alpha else ('exc' if np.median(d) > 0 or (np.median(d) == 0 and d.mean() > 0) else 'inh')
    mp, ms = pre.mean(), stim.mean()
    return dict(p=p, resp=resp, pre=mp, stim=ms, mi=(ms - mp) / (ms + mp) if ms + mp > 0 else 0.0,
                ratio=ms / mp if mp > 0 else np.nan)


def unit_response(spike_times, onsets, pre_s, on_s, alpha=ALPHA):
    pre = window_rates(spike_times, onsets, -pre_s, 0.0)
    stim = window_rates(spike_times, onsets, 0.0, on_s)
    return classify(pre, stim, alpha)


def profile(spike_times, onsets, t_from, t_to, bin_s):
    """Trial-averaged rate (Hz) in bins of bin_s from onset + t_from to onset + t_to; returns (bin centres, rate)."""
    edges = np.arange(t_from, t_to + 1e-9, bin_s)
    h = np.zeros(len(edges) - 1)
    st = np.asarray(spike_times)
    for o in onsets:
        h += np.histogram(st - o, bins=edges)[0]
    return (edges[:-1] + edges[1:]) / 2, h / len(onsets) / bin_s


def summarize(labels, ratios):
    labels, ratios = np.asarray(labels), np.asarray(ratios, float)
    out = {'n': int(len(labels))}
    for r in ('exc', 'inh', 'ns'):
        out[f'n_{r}'] = int((labels == r).sum())
        out[f'frac_{r}'] = float((labels == r).mean()) if len(labels) else np.nan
        sel = (labels == r) & np.isfinite(ratios)
        out[f'median_ratio_{r}'] = float(np.median(ratios[sel])) if sel.any() else np.nan
    return out
