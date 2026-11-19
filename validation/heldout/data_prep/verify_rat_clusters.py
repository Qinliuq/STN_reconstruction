"""
Check the shared metrics against the reference rat analysis before they are used on the model:
1. clustering port: unit_cluster_SSonly.csv features -> Ward k=2 labels (vs the saved ward_cluster), silhouette,
   Gaussian-null p95 / p (reference: 0.367 / 0.316 / ~0.002), and the PCA-2 GMM BIC (reference: k1 556, k2 583... see README)
2. feature port: features recomputed from the extracted spike trains vs the saved ones
3. sensitivity: the same tests after removing the duplicated unit copies (one row per file/label)
Writes results/heldout/H1/rat_verification.json and rat_features_recomputed.csv.

    python -m validation.heldout.data_prep.verify_rat_clusters
"""
import json
import os

import numpy as np
from sklearn.metrics import adjusted_rand_score

from validation.heldout import data_io
from validation.heldout.metrics import clustering, spike_features

OUT = os.path.join('results', 'heldout', 'H1')


def cluster_report(df):
    X = clustering.standardize(df)
    r = clustering.silhouette_null_test(X)
    r.update(clustering.bic_k1_k2(X))
    return r


def main():
    os.makedirs(OUT, exist_ok=True)
    ss = data_io.rat_ss_table()
    res = {}

    r = cluster_report(ss)
    labels = r.pop('labels')
    ref = ss.ward_cluster.values
    agree = max(np.mean(labels == ref), np.mean(labels == 1 - ref))
    res['clustering_73'] = dict(r, label_agreement=float(agree), ari=float(adjusted_rand_score(ref, labels)),
                                reference=dict(silhouette=0.367, null_p95=0.316, p=0.002, sizes=[33, 40], bic_k1=583, bic_k2=556, note='user-reported 556 vs 583 with k=2 preferred'))

    trains = data_io.rat_ss_trains()
    feats = spike_features.features_table(trains)
    feats[['file', 'label']] = feats.unit.tolist()
    feats = feats.drop(columns='unit')
    feats.to_csv(os.path.join(OUT, 'rat_features_recomputed.csv'), index=False)
    dedup = data_io.rat_dedup(ss)
    m = dedup.merge(feats, on=['file', 'label'], suffixes=('_ref', '_new'))
    diffs = {f: float(np.max(np.abs(m[f + '_ref'] - m[f + '_new']))) for f in spike_features.FEATURES + ['n_spikes']}
    res['feature_port'] = dict(n_units_compared=len(m), max_abs_diff=diffs)

    r = cluster_report(dedup)
    lab = r.pop('labels')
    res['clustering_dedup'] = dict(r, n=len(dedup), agreement_with_saved=float(max(np.mean(lab == dedup.ward_cluster.values),
                                                                                    np.mean(lab == 1 - dedup.ward_cluster.values))),
                                   saved_cluster_sizes=np.bincount(dedup.ward_cluster.values).tolist())
    with open(os.path.join(OUT, 'rat_verification.json'), 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))


if __name__ == '__main__':
    main()
