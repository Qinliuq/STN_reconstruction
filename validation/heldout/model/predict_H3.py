"""
H3 predictions: phase locking of model STN cells, before any comparison data is read.
- slow (0.5-2 Hz): spikes vs the phase of the cortical input (population arrival rate, band-passed, Hilbert; phase 0
  = peak of cortical input), in the parkinsonian + slow-oscillation runs (specs/H3_slowosc.py). Control: the same in
  the parkinsonian runs without slow oscillation (specs/H1_pd.py).
- beta (13-30 Hz), in the parkinsonian runs without slow oscillation (the analogue of activated-cortex epochs):
  (a) phase of the pallidal beta drive (the model's beta source; phase 0 = peak GPe rate), all time and within GPe
  beta bursts only; (b) the cortical input rate band-passed 13-30 Hz (the EEG analogue).
Cells are labelled by type and by their H1 cluster (H1/model_features_pd.csv, same seed and cell), so run
predict_H1_H2 first. Beta locking is also computed in the slow-oscillation runs (the analogue of slow-wave epochs).

    python -m validation.heldout.model.predict_H3
"""
import os

import numpy as np
import pandas as pd

from validation.heldout.metrics import phase
from validation.heldout.model import common

FS = 1000.0


def per_cell(runs, refs):
    rows = []
    for label, seed, d, _, _ in runs:
        sigs = refs(d)
        for i, typ, axis, st in common.cell_trains(d):
            row = dict(condition=label, seed=seed, cell=i, type=typ, axis=axis, n_spikes=len(st))
            for name, (ph, mask) in sigs.items():
                th = phase.spike_phases(st, ph, FS, 0.0, mask)
                ok = len(th) >= 50
                row.update({f'{name}_n': len(th), f'{name}_ppc': phase.ppc(th) if ok else np.nan,
                            f'{name}_pref': phase.preferred_phase(th) if ok else np.nan,
                            f'{name}_rayleigh_p': phase.rayleigh_p(th) if ok else np.nan})
            rows.append(row)
    return pd.DataFrame(rows)


def keep_after_t0(n):
    m = np.ones(n, bool)
    m[:int(common.T0_S * FS)] = False
    return m


def slow_refs(d):
    ctx = d['ctx_rate'].astype(float)
    return {'slow_ctx': (phase.band_phase(ctx, FS, 'slow'), keep_after_t0(len(ctx)))}


def beta_refs(d):
    ctx = d['ctx_rate'].astype(float)
    gph = d['gpe_beta_phase'].astype(float)
    m = keep_after_t0(len(gph))
    return {'beta_gpe_all': (gph, m), 'beta_gpe_burst': (gph, m & d['gpe_in_burst'].astype(bool)),
            'beta_ctx': (phase.band_phase(ctx, FS, 'beta'), keep_after_t0(len(ctx))),
            'slow_ctx_control': (phase.band_phase(ctx, FS, 'slow'), keep_after_t0(len(ctx)))}


def summarize(df, cols, by):
    out = {}
    for key, g in df.groupby(by):
        out[str(key)] = {c: dict(seedwise_median_ppc=common.mean_sd(g.groupby('seed')[f'{c}_ppc'].median()),
                                 frac_rayleigh_p01=common.mean_sd(g.groupby('seed')[f'{c}_rayleigh_p'].apply(lambda x: np.mean(x < 0.01))),
                                 circ_mean_pref=float(np.angle(np.nanmean(np.exp(1j * g[f'{c}_pref'].values)))))
                         for c in cols}
    return out


def main():
    out = common.out_dir('H3')
    so = common.load_runs('H3_slowosc')
    pd_runs = common.load_runs('H1_pd')
    clus = pd.read_csv(os.path.join(common.ROOT, 'results', 'heldout', 'H1', 'model_features_pd.csv'))[['seed', 'cell', 'bursty_seedwise']]
    so_refs = lambda d: {**slow_refs(d), **{k: v for k, v in beta_refs(d).items() if k != 'slow_ctx_control'}}
    a = per_cell(so, so_refs).merge(clus, on=['seed', 'cell'], how='left')
    b = per_cell(pd_runs, beta_refs).merge(clus, on=['seed', 'cell'], how='left')
    fa, fb = os.path.join(out, 'model_phase_slowosc.csv'), os.path.join(out, 'model_phase_pd_beta.csv')
    a.to_csv(fa, index=False)
    b.to_csv(fb, index=False)
    for df in (a, b):
        df['cluster'] = np.where(df.bursty_seedwise == True, 'bursty', np.where(df.bursty_seedwise == False, 'regular', 'excluded'))
    metrics = dict(
        slow=dict(by_type=summarize(a, ['slow_ctx'], 'type'), by_cluster=summarize(a, ['slow_ctx'], 'cluster'),
                  all=summarize(a.assign(all='all'), ['slow_ctx'], 'all')),
        beta_in_slowosc_runs=summarize(a.assign(all='all'), ['beta_ctx', 'beta_gpe_all', 'beta_gpe_burst'], 'all'),
        beta=dict(by_type=summarize(b, ['beta_gpe_all', 'beta_gpe_burst', 'beta_ctx', 'slow_ctx_control'], 'type'),
                  by_cluster=summarize(b, ['beta_gpe_all', 'beta_gpe_burst', 'beta_ctx', 'slow_ctx_control'], 'cluster'),
                  all=summarize(b.assign(all='all'), ['beta_gpe_all', 'beta_gpe_burst', 'beta_ctx', 'slow_ctx_control'], 'all')),
        phase_convention='phase 0 = peak of the reference (cortical input rate or GPe rate); radians in (-pi, pi]')
    common.write_prediction('H3', metrics, so + pd_runs, [fa, fb])


if __name__ == '__main__':
    main()
