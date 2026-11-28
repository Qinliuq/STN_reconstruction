"""
Loaders for the held-out experimental data. Only data_prep/ and compare/ may import this module: the model
runs and predictions (model/) must not see the data (enforced by tests/test_isolation.py).

The data stay outside the repository; set HELDOUT_DATA_DIR to the unpacked heldout_bundle folder.
"""
import glob
import os

import numpy as np
import pandas as pd

DATA_DIR = os.environ.get('HELDOUT_DATA_DIR', '/home/user/heldout_data/heldout_bundle')
DUP_TOL_S = 1e-6 # spike times closer than this within one label are the same spike (see rat_ss_trains)


def path(*p):
    return os.path.join(DATA_DIR, *p)


# ------------------------------------------------ rat 6-OHDA ------------------------------------------------
def rat_ss_table():
    """unit_cluster_SSonly.csv: the 73 SS rows with the saved features, ward_cluster (0 bursty, 1 regular) and PC1."""
    return pd.read_csv(path('rat_6ohda', 'unit_cluster_SSonly.csv'))


def rat_ss_trains():
    """{(file, label): spike times (s)} for the SS units of the extractor output (rat_mat/spikes_*.csv).

    The Spike2 exports contain 28 (file, label) titles twice: the second copy is the same unit (identical, or a
    subset of the first copy's spike times). The extractor writes both copies under one label, so here coincident
    spikes (< DUP_TOL_S apart) are merged, which recovers the larger copy."""
    out = {}
    for f in sorted(glob.glob(path('rat_mat', 'spikes_*.csv'))):
        sp = pd.read_csv(f)
        sp = sp[sp.unit_type == 'SS']
        for (fname, label), g in sp.groupby(['file', 'label']):
            t = np.sort(g.spike_time_s.values)
            out[(fname, label)] = t[np.r_[True, np.diff(t) > DUP_TOL_S]]
    return out


def rat_eeg(session):
    """(signal, fs, t0_s, channel) for a session ('NPR040d'): IpsiEEG, falling back to ContraEEG; None if absent."""
    d = np.load(path('rat_mat', f'eeg_{session}.npz'))
    for ch in ('IpsiEEG', 'ContraEEG'):
        if ch in d.files:
            return d[ch].astype(np.float64), float(d['fs']), float(d['t0_s']), ch
    return None


def rat_dedup(df):
    """Keep one row per (file, label): the copy with more spikes (the others are copies of the same unit)."""
    return df.sort_values('n_spikes', ascending=False).drop_duplicates(['file', 'label']).sort_index()


# ------------------------------------------------ mouse ChR2 ------------------------------------------------
def chr2_binned():
    return pd.read_csv(path('chr2_mouse', 'analysis_stn_chr2', 'binned_response.csv'))


def chr2_d2_wc():
    return pd.read_csv(path('chr2_mouse', 'analysis_stn_chr2', 'd2_wc_response_xtalk_removed.csv'))
