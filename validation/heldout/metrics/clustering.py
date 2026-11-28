"""
Two-cluster test on the 8 spike-train features, shared by data and model.

Port of rat_6ohda/cluster_test.py: StandardScaler fitted within the analysed subset, Ward agglomerative k = 2,
silhouette; null = 1000 draws from a multivariate Gaussian with the subset's covariance (of the standardized
features), clustered the same way, rng seed 0; p = fraction of null silhouettes >= observed.
BIC check (not in the original scripts, reimplemented from the description): GaussianMixture on the first 2 PCs
of the standardized features, k = 1 vs k = 2, full covariance, several inits; lower BIC is better.
"""
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

from .spike_features import FEATURES


def standardize(df, cols=FEATURES):
    return StandardScaler().fit_transform(df[cols].values.astype(float))


def ward2(X):
    return AgglomerativeClustering(n_clusters=2, linkage='ward').fit_predict(X)


def silhouette_null_test(X, n_null=1000, seed=0):
    n, p = X.shape
    labels = ward2(X)
    sil = silhouette_score(X, labels)
    cov = np.cov(X, rowvar=False)
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_null):
        Xn = rng.multivariate_normal(np.zeros(p), cov, size=n)
        ln = ward2(Xn)
        if len(set(ln)) == 2:
            null.append(silhouette_score(Xn, ln))
    null = np.array(null)
    return dict(labels=labels, silhouette=float(sil), null_mean=float(null.mean()),
                null_p95=float(np.percentile(null, 95)), p=float(np.mean(null >= sil)),
                sizes=np.bincount(labels).tolist())


def bic_k1_k2(X, n_pcs=2, n_init=10, seed=0):
    Z = PCA(n_components=n_pcs).fit_transform(X)
    out = {}
    for k in (1, 2):
        gm = GaussianMixture(n_components=k, covariance_type='full', n_init=n_init, random_state=seed).fit(Z)
        out[f'bic_k{k}'] = float(gm.bic(Z))
    out['prefers_k2'] = out['bic_k2'] < out['bic_k1']
    return out


def bursty_cluster(df, labels, feature='cv'):
    """Index (0/1) of the cluster with the higher median of `feature` (the bursty/pausy cluster)."""
    import pandas as pd
    med = pd.Series(df[feature].values).groupby(labels).median()
    return int(med.idxmax())


def assign_to_reference(X_new, X_ref, ref_labels):
    """Nearest-centroid assignment of new (already standardized in the reference space) points to reference clusters."""
    cents = np.array([X_ref[ref_labels == k].mean(0) for k in (0, 1)])
    d = ((X_new[:, None, :] - cents[None]) ** 2).sum(-1)
    return d.argmin(1)
