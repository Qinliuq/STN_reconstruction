"""Unit tests of the shared metrics on synthetic spike trains with known answers.

    python -m pytest validation/heldout/tests -q
"""
import os

import numpy as np
import pandas as pd
import pytest
from scipy.special import i0, i1

from validation.heldout.metrics import clustering, phase, response, spike_features as sf

RNG = np.random.default_rng(0)


def poisson(rate, T, rng=RNG):
    return np.sort(rng.uniform(0, T, rng.poisson(rate * T)))


def regular(isi, T, jitter=0.01, rng=RNG):
    t = np.arange(0, T, isi)
    return np.sort(t + rng.normal(0, jitter * isi, len(t)))


def bursty(T, burst_rate=2.0, n_in=5, intra=0.005, rng=RNG):
    starts = poisson(burst_rate, T, rng)
    starts = starts[np.r_[True, np.diff(starts) > n_in * intra * 3]]
    return np.sort(np.concatenate([s + np.arange(n_in) * intra for s in starts]))


def test_poisson_features():
    f = sf.unit_features(poisson(20, 300))
    assert abs(f['rate_hz'] - 20) < 1
    assert abs(f['cv'] - 1) < 0.05 and abs(f['lv'] - 1) < 0.05
    assert abs(f['short_isi_rate_relative'] - (1 - np.exp(-0.5))) < 0.02
    assert abs(f['pause_time_frac'] - 5 * np.exp(-4)) < 0.015
    assert abs(f['frac_isi_lt10ms'] - (1 - np.exp(-0.2))) < 0.02
    # S >= 3 has false positives over long Poisson trains: ~0.10 of spikes at 20 Hz / 300 s (Monte Carlo), so a
    # burst share well below 0.1 (rat 'regular' cluster: 0.0015) means sub-Poisson, more regular firing
    assert 0.05 < f['burst_share'] < 0.15
    assert f['log_isi_bimodality'] < 0.555 # unimodal (log-ISI of a Poisson train is skewed but unimodal)


def test_regular_features():
    f = sf.unit_features(regular(0.05, 100))
    assert f['cv'] < 0.03 and f['lv'] < 0.01
    assert f['short_isi_rate_relative'] == 0 and f['pause_time_frac'] == 0 and f['burst_share'] == 0


def test_bursty_features():
    f = sf.unit_features(bursty(300))
    assert f['burst_share'] > 0.8
    assert f['cv'] > 1.5 and f['pause_time_frac'] > 0.6


def test_bimodality_coefficient():
    rng = np.random.default_rng(4)
    # the reference code puts the Pearson kurtosis (normal = 3) in the denominator of Sarle's coefficient, where
    # Sarle uses the excess kurtosis; values are ~half the usual ones: normal 1/6 (usual 1/3), uniform 1/4.8
    # (usual 5/9), a symmetric two-point mixture -> 1/4 (usual 1). Kept as in the reference.
    assert abs(sf.bimodality_coeff(rng.normal(size=20000)) - 1 / 6) < 0.01
    assert abs(sf.bimodality_coeff(rng.uniform(size=20000)) - 1 / 4.8) < 0.01
    assert sf.bimodality_coeff(np.r_[rng.normal(-3, 0.5, 5000), rng.normal(3, 0.5, 5000)]) > 0.23


def test_inclusion_criteria():
    assert sf.unit_features(np.arange(99) * 0.1) is None
    assert sf.unit_features(np.repeat(np.arange(60) * 0.1, 2)) is None # 120 spikes but only 59 positive ISIs
    assert sf.unit_features(np.arange(100) * 0.1) is not None


def test_surprise_detects_planted_burst():
    t = regular(0.1, 50, jitter=0.0)
    t = np.sort(np.r_[t, 25.0 + 0.001 + np.arange(8) * 0.004])
    isi = np.diff(t)
    mask = sf.poisson_surprise_bursts(t, isi, len(t) / (t[-1] - t[0]))
    assert mask[(t > 25.0) & (t < 25.04)].all() and mask.sum() <= 11


def test_clustering_two_groups_vs_one():
    rng = np.random.default_rng(1)
    feats = sf.FEATURES
    two = pd.DataFrame(np.r_[rng.normal(0, 1, (40, 8)), rng.normal(4, 1, (40, 8))], columns=feats)
    one = pd.DataFrame(rng.normal(0, 1, (80, 8)), columns=feats)
    r2 = clustering.silhouette_null_test(clustering.standardize(two), n_null=200)
    r1 = clustering.silhouette_null_test(clustering.standardize(one), n_null=200)
    assert r2['p'] < 0.01 and sorted(r2['sizes']) == [40, 40]
    assert r1['p'] > 0.05
    assert clustering.bic_k1_k2(clustering.standardize(two))['prefers_k2']
    assert not clustering.bic_k1_k2(clustering.standardize(one))['prefers_k2']


@pytest.mark.parametrize('kappa,mu', [(1.0, 0.5), (2.0, -2.0)])
def test_phase_locking_von_mises(kappa, mu):
    fs, T, f0 = 1000.0, 400.0, 1.0
    t = np.arange(int(T * fs)) / fs
    ref = np.cos(2 * np.pi * f0 * t)
    rng = np.random.default_rng(2)
    cycles = np.arange(5, int(T * f0) - 5)
    th = rng.vonmises(mu, kappa, size=len(cycles) * 5)
    spikes = np.repeat(cycles, 5) / f0 + np.mod(th, 2 * np.pi) / (2 * np.pi * f0)
    r = phase.locking(spikes, ref, fs, 'slow')
    expected = (i1(kappa) / i0(kappa)) ** 2
    assert abs(r['ppc'] - expected) < 0.05
    assert abs(np.angle(np.exp(1j * (r['pref_phase'] - mu)))) < 0.15
    assert r['rayleigh_p'] < 1e-6


def test_phase_locking_uniform():
    fs = 1000.0
    t = np.arange(int(200 * fs)) / fs
    r = phase.locking(poisson(10, 200), np.cos(2 * np.pi * 20 * t), fs, 'beta')
    assert abs(r['ppc']) < 0.005


def test_response_classification():
    rng = np.random.default_rng(3)
    onsets = 5.0 + np.arange(10) * 5.0
    pieces = []
    for o in onsets: # 20 Hz pre / post, 10 Hz during 2 s on
        pieces += [o - 1 + np.sort(rng.uniform(0, 1, rng.poisson(20))), o + np.sort(rng.uniform(0, 2, rng.poisson(20))),
                   o + 2 + np.sort(rng.uniform(0, 2, rng.poisson(40)))]
    st = np.concatenate(pieces)
    r = response.unit_response(st, onsets, 1.0, 2.0)
    assert r['resp'] == 'inh' and abs(r['ratio'] - 0.5) < 0.15
    r0 = response.unit_response(poisson(20, 60, rng), onsets, 1.0, 2.0)
    assert r0['resp'] == 'ns'
    c, p = response.profile(st, onsets, -1, 4, 1.0)
    assert len(c) == 5 and p[1] < p[0]


def test_burst_baseline():
    from validation.heldout.metrics.burst_baseline import burst_pause
    rng = np.random.default_rng(3)
    T = 10.0
    # Poisson trains pass only at about the 5% level of the rule (burst share > p95 and CV > 1: ~2-3%)
    fp = np.mean([burst_pause(np.sort(rng.uniform(0, T, 200)), T, seed=i)['burst_pause'] for i in range(40)])
    assert fp <= 0.1
    regular = np.arange(0.01, T, 0.05)
    assert not burst_pause(regular, T)['burst_pause']
    # 3-spike bursts (5 ms ISI) at ~5 Hz with jittered inter-burst intervals
    starts = np.cumsum(rng.uniform(0.15, 0.25, 50))
    bursty = np.sort(np.concatenate([starts + k * 0.005 for k in range(3)]))
    bursty = bursty[bursty < T]
    r = burst_pause(bursty, T)
    assert r['burst_pause'] and r['burst_share'] > r['p95_poisson'] and r['cv'] > 1
    assert not burst_pause(regular[:20], T)['burst_pause'] # too few spikes
