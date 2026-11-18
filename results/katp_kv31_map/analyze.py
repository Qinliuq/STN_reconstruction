"""
Kv3.1 x K_ATP single-cell map (PV+): per-cell measures and admissibility (criteria.md, fixed before the runs).

    python results/katp_kv31_map/analyze.py

Reads raw/kv<Kv3.1>_kir<Kir>.npz (11 s default protocols) and raw/..._depol.npz (4 s depolarizing steps), written by
src/single_cell.py --save. Writes map.csv, diagnostic.csv and map.png.
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'src'))
from validation.heldout.metrics.burst_baseline import burst_pause  # noqa: E402

T0, MS = 1000.0, 1e-3
SPONT, HOLD, REB = 'spontaneous (0 nA)', 'hold -0.16 nA', 'rebound (-0.25 nA, 1000-1500 ms)'
DEPOL = 'depolarizing +0.8 nA'


def find_bursts(spikes, max_isi=20.0, min_spikes=3, pause_ratio=3.0):
    # same definition as src/single_cell.find_bursts (copied: single_cell.py runs a simulation on import)
    isi = np.diff(spikes)
    runs, run = [], [0]
    for i, x in enumerate(isi):
        if x < max_isi:
            run.append(i + 1)
        else:
            runs.append(run)
            run = [i + 1]
    runs.append(run)
    out = []
    for run in runs:
        if len(run) < min_spikes or run[0] == 0 or run[-1] == len(spikes) - 1:
            continue
        intra = np.mean(np.diff(spikes[run]))
        if min(isi[run[0] - 1], isi[run[-1]]) > pause_ratio * intra:
            out.append(spikes[run])
    return out


def window_stats(sp, t0, t1):
    w = sp[(sp >= t0) & (sp < t1)]
    isi = np.diff(w)
    b = [x for x in find_bursts(sp) if t0 <= x[0] < t1]
    return dict(rate=1000 * len(w) / (t1 - t0), cv=isi.std() / isi.mean() if len(isi) >= 2 else np.nan,
                bursts_per_s=1000 * len(b) / (t1 - t0), spikes_per_burst=np.mean([len(x) for x in b]) if b else np.nan)


def rebound(sp, release=1500.0, max_isi=30.0):
    post = sp[sp >= release]
    if len(post) < 2:
        return len(post), 0.0
    isi = np.diff(post)
    end = int(np.argmax(isi > max_isi)) if (isi > max_isi).any() else len(isi)
    n, dur = end + 1, post[end] - post[0]
    return n, 1000 * (n - 1) / dur if dur > 0 else 0.0


def load(path):
    d = np.load(path, allow_pickle=False)
    labels, prots = list(d['labels']), list(d['protocols'])
    gbars = [float(re.search(r'katp_gbar=([0-9.e-]+)', l).group(1)) for l in labels]
    return d, gbars, prots


rows = []
for f in sorted(glob.glob(os.path.join(HERE, 'raw', 'kv*_kir*.npz'))):
    if f.endswith('_depol.npz'):
        continue
    kv, kir = map(float, re.search(r'kv([0-9.]+)_kir([0-9.]+)\.npz', f).groups())
    d, gbars, prots = load(f)
    dd, _, dprots = load(f.replace('.npz', '_depol.npz'))
    dur = float(d['duration'])
    for c, g in enumerate(gbars):
        sp = lambda p: d[f'spk_{c}_{prots.index(p)}']
        s = sp(SPONT)
        early = window_stats(s, 1000, 4000)  # window of the 2026-10-18 run
        full = window_stats(s, T0, dur)
        bp = burst_pause(s[s >= T0] * MS, (dur - T0) * MS, seed=c)
        h = window_stats(sp(HOLD), T0, dur)
        rn, rf = rebound(sp(REB))
        i8 = dprots.index(DEPOL)
        s8 = dd[f'spk_{c}_{i8}']
        dDur = float(dd['duration'])
        last8 = s8[-1] if len(s8) else 0.0
        trough8 = float(dd[f'vminlate_{c}_{i8}'])
        rows.append(dict(
            kv31=kv, kir=kir, katp_gbar=g,
            spont_rate=full['rate'], spont_cv=full['cv'], spont_bursts_per_s=full['bursts_per_s'],
            spont_spikes_per_burst=full['spikes_per_burst'],
            spont_rate_1_4s=early['rate'], spont_cv_1_4s=early['cv'], spont_bursts_per_s_1_4s=early['bursts_per_s'],
            spont_spikes_per_burst_1_4s=early['spikes_per_burst'],
            burst_share=bp['burst_share'], p95_poisson=bp['p95_poisson'], n_spikes=bp['n_spikes'],
            hold016_rate=h['rate'], hold016_bursts_per_s=h['bursts_per_s'], hold016_spikes_per_burst=h['spikes_per_burst'],
            rebound_spikes=rn, rebound_hz=rf,
            depol08_rate=1000 * len(s8[s8 >= T0]) / (dDur - T0), depol08_last_spike=last8, depol08_trough_late=trough8,
            c1_burst_pause=bp['burst_pause'],
            c2_no_block=bool(last8 >= dDur - 500 and trough8 < -50),
            c3_phenotype=bool(h['bursts_per_s'] >= 1 and (h['spikes_per_burst'] or 0) >= 3 and rn >= 8 and rf >= 100),
        ))

df = pd.DataFrame(rows).sort_values(['kir', 'kv31', 'katp_gbar'])
df['admissible'] = df.c1_burst_pause & df.c2_no_block & df.c3_phenotype
diag = df[((df.kv31 == 1.0) | (df.kv31 == 1.2)) & df.kir.isin([0.3, 0.5])]
diag.to_csv(os.path.join(HERE, 'diagnostic.csv'), index=False, float_format='%.4g')
grid = df[df.kir == 0.3]
grid.to_csv(os.path.join(HERE, 'map.csv'), index=False, float_format='%.4g')

# ---- figure: map over Kv3.1 x gbar (Kir 0.3)
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

kvs = sorted(grid.kv31.unique())
gs = sorted(grid.katp_gbar.unique())
fig, axs = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
for ax, (col, title, cmap) in zip(axs, [('spont_cv', 'spontaneous ISI CV', 'viridis'),
                                        ('spont_rate', 'spontaneous rate (Hz)', 'magma'),
                                        ('burst_share', 'Poisson-surprise burst share', 'cividis')]):
    M = np.array([[grid[(grid.kv31 == k) & (grid.katp_gbar == g)][col].iloc[0] for g in gs] for k in kvs])
    im = ax.imshow(M, origin='lower', aspect='auto', cmap=cmap)
    fig.colorbar(im, ax=ax, shrink=0.8)
    for i, k in enumerate(kvs):
        for j, g in enumerate(gs):
            r = grid[(grid.kv31 == k) & (grid.katp_gbar == g)].iloc[0]
            flags = ('B' if r.c1_burst_pause else '-') + ('D' if r.c2_no_block else 'x') + ('P' if r.c3_phenotype else '-')
            ax.text(j, i, f'{r[col]:.2g}\n{flags}', ha='center', va='center', fontsize=7,
                    color='white', fontweight='bold' if r.admissible else 'normal')
            if r.admissible:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec='red', lw=2.5))
    ax.set_xticks(range(len(gs)), [f'{g * 1e5:g}' for g in gs])
    ax.set_yticks(range(len(kvs)), [f'{k:g}' for k in kvs])
    ax.set_xlabel('K_ATP gbar (x1e-5 S/cm2)')
    ax.set_ylabel('PV+ Kv3.1 scale (Kir 0.3)')
    ax.set_title(title, fontsize=10)
fig.suptitle('PV+ Kv3.1 x K_ATP (k_Ca 20 uM). Flags: B spontaneous burst/pause, D no block at +0.8 nA, '
             'P -0.16 nA bursts + rebound kept; red box = admissible', fontsize=9)
fig.savefig(os.path.join(HERE, 'map.png'), dpi=130)

pd.set_option('display.width', 250)
cols = ['kv31', 'kir', 'katp_gbar', 'spont_rate', 'spont_cv', 'spont_spikes_per_burst', 'burst_share', 'p95_poisson',
        'hold016_bursts_per_s', 'hold016_spikes_per_burst', 'rebound_spikes', 'rebound_hz', 'depol08_rate',
        'depol08_trough_late', 'c1_burst_pause', 'c2_no_block', 'c3_phenotype', 'admissible']
print(df[cols].to_string(index=False, float_format=lambda x: f'{x:.3g}'))
