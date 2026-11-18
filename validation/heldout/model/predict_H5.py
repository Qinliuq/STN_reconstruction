"""
H5 prediction: healthy -> parkinsonian change in configuration D, written before the primate data are downloaded.

Healthy: specs/H5_healthy.py (5 seeds x 20.5 s). Parkinsonian: the H1 runs (specs/H1_pd.py), checked against the spike
hashes recorded in results/heldout/H1/prediction.json. Per cell and seed: metrics.unit_summary.h5_summary over
t >= 0.5 s (rate, CV, LV, pause-time fraction, frac ISI < 10 ms, burst share, oscillation test, peak frequency).
Per seed: population medians per state, the PD / healthy rate ratio and PD - healthy differences of CV, LV and
pause-time fraction, oscillatory fractions and the share of PD oscillatory cells peaking at 4-17 Hz. All also by
PV+ / PV-.

    python -m validation.heldout.model.predict_H5
"""
import json
import os

import numpy as np
import pandas as pd

from validation.heldout.metrics.unit_summary import h5_summary, FEATURES
from validation.heldout.model import common

T0 = common.T0_S


def check_h1_hashes():
    p = json.load(open(os.path.join(common.ROOT, 'results', 'heldout', 'H1', 'prediction.json'), encoding='utf-8'))
    for r in p['runs']:
        if r['label'] != 'pd':
            continue
        f = os.path.join(common.ROOT, *r['spikes'].replace('\\', '/').split('/'))
        if common.sha256(f) != r['spikes_sha256']:
            raise SystemExit(f'{f} differs from the H1 record; re-run the PD runs')


def cells(runs, state):
    rows, excluded = [], []
    for label, seed, d, cfg, path in runs:
        for i, typ, axis, st in common.cell_trains(d):
            s = h5_summary(st, T0, T_END, seed=1000 * seed + i)
            if s is None:
                excluded.append(dict(state=state, seed=seed, cell=i, type=typ, n_spikes=len(st)))
                continue
            rows.append(dict(state=state, seed=seed, cell=i, type=typ, axis=axis, **s))
    return rows, excluded


def per_seed(df):
    out = []
    for grp, sub in [('all', df), ('PV+', df[df.type == 'PV+']), ('PV-', df[df.type == 'PV-'])]:
        for seed in sorted(df.seed.unique()):
            h = sub[(sub.state == 'healthy') & (sub.seed == seed)]
            p = sub[(sub.state == 'parkinsonian') & (sub.seed == seed)]
            r = dict(group=grp, seed=seed, n_healthy=len(h), n_pd=len(p))
            for k in FEATURES:
                r[f'{k}_healthy'] = h[k].median()
                r[f'{k}_pd'] = p[k].median()
            r['rate_ratio'] = r['rate_hz_pd'] / r['rate_hz_healthy']
            for k in ('cv', 'lv', 'pause_time_frac'):
                r[f'{k}_diff'] = r[f'{k}_pd'] - r[f'{k}_healthy']
            r['osc_frac_healthy'] = h.oscillatory.mean()
            r['osc_frac_pd'] = p.oscillatory.mean()
            po = p[p.oscillatory]
            r['pd_peak_4_17_frac'] = float(((po.peak_hz >= 4) & (po.peak_hz < 17)).mean()) if len(po) else np.nan
            r['pd_peak_median_hz'] = po.peak_hz.median() if len(po) else np.nan
            out.append(r)
    return pd.DataFrame(out)


def main():
    global T_END
    spec_ms = common.PRIORS['H5_primate_change']['duration_ms']
    T_END = spec_ms / 1000
    check_h1_hashes()
    runs_h = common.load_runs('H5_healthy')
    runs_p = [r for r in common.load_runs('H1_pd') if r[0] == 'pd']
    rh, eh = cells(runs_h, 'healthy')
    rp, ep = cells(runs_p, 'parkinsonian')
    df = pd.DataFrame(rh + rp)
    out = common.out_dir('H5')
    f_cells = os.path.join(out, 'model_cells.csv')
    df.to_csv(f_cells, index=False, float_format='%.6g')
    ps = per_seed(df)
    f_seed = os.path.join(out, 'model_per_seed.csv')
    ps.to_csv(f_seed, index=False, float_format='%.6g')
    metrics = dict(n_cells=dict(healthy=len(rh), parkinsonian=len(rp)), excluded=eh + ep, window_s=[T0, T_END])
    for grp in ('all', 'PV+', 'PV-'):
        g = ps[ps.group == grp]
        metrics[grp] = {k: common.mean_sd(g[k]) for k in
                        ['rate_ratio', 'cv_diff', 'lv_diff', 'pause_time_frac_diff', 'osc_frac_healthy', 'osc_frac_pd',
                         'pd_peak_4_17_frac', 'pd_peak_median_hz'] + [f'{k}_{s}' for k in FEATURES for s in ('healthy', 'pd')]}
    common.write_prediction('H5', metrics, runs_h + runs_p, files=[f_cells, f_seed])
    a = metrics['all']
    print(f"rate ratio {a['rate_ratio']['mean']:.3f} +- {a['rate_ratio']['sd']:.3f}; osc frac healthy "
          f"{a['osc_frac_healthy']['mean']:.2f} -> PD {a['osc_frac_pd']['mean']:.2f}; PD peaks 4-17 Hz "
          f"{a['pd_peak_4_17_frac']['mean']:.2f}")


T_END = None

if __name__ == '__main__':
    main()
