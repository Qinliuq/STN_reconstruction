"""
H1 / H2 predictions from the parkinsonian runs (specs/H1_pd.py), before any comparison data is read.
H1: the 8 reference features per cell; two-cluster test (Ward + silhouette vs Gaussian null, PCA-2 GMM BIC) per seed
and pooled; bursty-cluster proportion; bursty cluster vs PV+ (Cramer's V, Fisher). Secondary (labelled): the same on
the parkinsonian + slow-oscillation runs (specs/H3_slowosc.py), if they are finished.
H2: pooled per-cell rate / CV / burst share (written to H2/prediction.json).

    python -m validation.heldout.model.predict_H1_H2
"""
import os

import numpy as np
import pandas as pd

from validation.heldout.metrics import clustering, spike_features as sf, stats
from validation.heldout.model import common


def features(runs):
    rows, excluded = [], []
    for label, seed, d, _, _ in runs:
        for i, typ, axis, st in common.cell_trains(d):
            f = sf.unit_features(st)
            if f is None:
                excluded.append(dict(seed=seed, cell=i, type=typ, n_spikes=len(st)))
            else:
                rows.append(dict(condition=label, seed=seed, cell=i, type=typ, axis=axis, **f))
    return pd.DataFrame(rows), excluded


def cluster_block(df):
    X = clustering.standardize(df)
    r = clustering.silhouette_null_test(X)
    labels = r.pop('labels')
    r.update(clustering.bic_k1_k2(X))
    b = clustering.bursty_cluster(df, labels)
    bursty = labels == b
    r['bursty_frac'] = float(bursty.mean())
    r['bursty_n'] = int(bursty.sum())
    cv = stats.cramers_v(np.where(bursty, 'bursty', 'regular'), df['type'].values)
    r['vs_type'] = cv
    r['frac_PVp_in_bursty'] = float((df['type'].values[bursty] == 'PV+').mean()) if bursty.any() else np.nan
    r['frac_bursty_in_PVp'] = float(bursty[df['type'].values == 'PV+'].mean())
    r['cluster_medians'] = {k: df.loc[bursty == (k == 'bursty'), sf.FEATURES].median().to_dict() for k in ('bursty', 'regular')}
    return r, bursty


def analyse(df):
    per_seed, flags = {}, np.zeros(len(df), bool)
    for seed, g in df.groupby('seed'):
        r, b = cluster_block(g)
        per_seed[int(seed)] = r
        flags[g.index.values] = b
    pooled, pooled_b = cluster_block(df)
    summ = {k: common.mean_sd([per_seed[s][k] for s in per_seed]) for k in
            ('silhouette', 'null_p95', 'p', 'bic_k1', 'bic_k2', 'bursty_frac', 'frac_PVp_in_bursty', 'frac_bursty_in_PVp')}
    summ['V'] = common.mean_sd([per_seed[s]['vs_type'].get('V', np.nan) for s in per_seed])
    summ['n_seeds_sig_p05'] = int(sum(per_seed[s]['p'] < 0.05 for s in per_seed))
    summ['n_seeds_bic_k2'] = int(sum(per_seed[s]['prefers_k2'] for s in per_seed))
    return dict(per_seed=per_seed, pooled=pooled, summary=summ), flags, pooled_b


def main():
    out = common.out_dir('H1')
    runs = common.load_runs('H1_pd')
    df, excl = features(runs)
    res, flags, pooled_b = analyse(df)
    df['bursty_seedwise'] = flags
    df['bursty_pooled'] = pooled_b
    fpath = os.path.join(out, 'model_features_pd.csv')
    df.to_csv(fpath, index=False)
    metrics = dict(primary_pd=res, excluded_cells=excl, n_cells=int(len(df)),
                   pv_fraction=float((df.type == 'PV+').mean()))
    files = [fpath]
    allruns = list(runs)
    try:
        so = common.load_runs('H3_slowosc')
        df2, excl2 = features(so)
        res2, f2, p2 = analyse(df2)
        df2['bursty_seedwise'], df2['bursty_pooled'] = f2, p2
        f2path = os.path.join(out, 'model_features_pd_slowosc.csv')
        df2.to_csv(f2path, index=False)
        metrics['secondary_pd_slowosc'] = dict(res2, excluded_cells=excl2)
        files.append(f2path)
        allruns += so
    except SystemExit as e:
        print('secondary (slow oscillation) skipped:', e)
    common.write_prediction('H1', metrics, allruns, files)

    # H2: pooled per-cell distributions (no fitting, sanity check)
    h2 = {k: dict(all=df[k].tolist(), **{t: df.loc[df.type == t, k].tolist() for t in ('PV+', 'PV-')},
                  median=float(df[k].median())) for k in ('rate_hz', 'cv', 'burst_share')}
    common.write_prediction('H2', dict(pooled_check=True, per_cell=h2), runs, [fpath])


if __name__ == '__main__':
    main()
