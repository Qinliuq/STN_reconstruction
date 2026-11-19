"""
H4 predictions: virtual D2- / D1-MSN stimulation in the healthy state (specs/H4_opto.py), before any data is read.
Per cell and seed: Wilcoxon of the 2 s stim rate vs the 1 s pre rate across the 10 cycles (p < 0.01) -> exc / inh / ns,
ratio stim / pre; window ratios (onset 0-0.25 s, late 0.25-2 s, post 0-2 s after offset; vs pre); trial-averaged
profiles in 0.25 s and 1 s bins; responder enrichment by PV type (Fisher).

    python -m validation.heldout.model.predict_H4
"""
import os

import numpy as np
import pandas as pd

from validation.heldout.metrics import response, stats
from validation.heldout.model import common

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src'))
from input_modulation import stim_onsets  # noqa: E402 (same schedule as the model's input hook)

P = common.PRIORS['H4_optogenetics']
PROTO = P['protocol']


def main():
    out = common.out_dir('H4')
    runs = common.load_runs('H4_opto')
    rows, prof = [], []
    for label, seed, d, _, _ in runs:
        on = stim_onsets(float(d['t1']) / 1000, **PROTO)
        on = on[:PROTO['n_cycles']]
        for i, typ, axis, st in common.cell_trains(d):
            r = response.unit_response(st, on, PROTO['pre_s'], PROTO['on_s'])
            pre = r['pre']
            w = lambda a, b: response.window_rates(st, on, a, b).mean()
            rows.append(dict(scenario=label, seed=seed, cell=i, type=typ, axis=axis, **r,
                             onset_ratio=w(0, 0.25) / pre if pre > 0 else np.nan,
                             late_ratio=w(0.25, PROTO['on_s']) / pre if pre > 0 else np.nan,
                             post_ratio=w(PROTO['on_s'], PROTO['on_s'] + PROTO['post_s']) / pre if pre > 0 else np.nan))
            for bw in (0.25, 1.0):
                c, p = response.profile(st, on, -PROTO['pre_s'], PROTO['on_s'] + PROTO['post_s'], bw)
                prof += [dict(scenario=label, seed=seed, cell=i, type=typ, bin_s=bw, t=float(ci), rate=float(pi), pre=pre)
                         for ci, pi in zip(c, p)]
    df = pd.DataFrame(rows)
    pf = pd.DataFrame(prof)
    f1, f2 = os.path.join(out, 'model_responses.csv'), os.path.join(out, 'model_profiles.csv')
    df.to_csv(f1, index=False)
    pf.to_csv(f2, index=False)
    metrics = {}
    for sc, g in df.groupby('scenario'):
        per_seed = {int(s): response.summarize(gs.resp.values, gs.ratio.values) for s, gs in g.groupby('seed')}
        m = {k: common.mean_sd([per_seed[s][k] for s in per_seed]) for k in
             ('frac_exc', 'frac_inh', 'median_ratio_exc', 'median_ratio_inh')}
        m['pop_ratio'] = common.mean_sd(g.groupby('seed').apply(lambda x: x.stim.sum() / x.pre.sum()))
        m['sign'] = 'excitation' if m['pop_ratio']['mean'] > 1 else 'inhibition'
        m['pooled'] = response.summarize(g.resp.values, g.ratio.values)
        for k in ('onset_ratio', 'late_ratio', 'post_ratio'):
            m[k] = common.mean_sd(g.groupby('seed')[k].median())
        resp = np.where(g.resp != 'ns', 'responder', 'ns')
        m['type_enrichment_pooled'] = stats.cramers_v(resp, g.type.values)
        m['responder_frac_by_type'] = {t: common.mean_sd(gt.groupby('seed').resp.apply(lambda x: np.mean(x != 'ns')))
                                       for t, gt in g.groupby('type')}
        m['label'] = P['scenarios'][sc]['label']
        metrics[sc] = dict(summary=m, per_seed=per_seed)
    common.write_prediction('H4', metrics, runs, [f1, f2])


if __name__ == '__main__':
    main()
